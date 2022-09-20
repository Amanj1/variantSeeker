#!/usr/bin/env nextflow

/*
How to run:
nextflow -C variantSeeker.config run variantSeeker.nf -profile amanj
*/

/* input files */
//fastq_files = Channel.fromFilePairs("${params.fastq_dir}/*.fq.gz"){ (it.name =~ /P[0-9]{3,5}_[0-9]{3,5}/)[0]}
fastq_files = Channel.fromFilePairs("${params.fastq_dir}/**/*_{1,2,unpaired}.fq.gz",size:3)
html_files = Channel.fromFilePairs("${params.html_dir}/html_{start,end}.txt",size:2)

fastq_files.into{reformat_PE_in;
reformat_S_in;
extraction_in}

html_files.into{html_each;
html_all}
/**
Nextflow processes
**/

process reformat_fastq_reads_PE{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reformation", mode:'link'

  input:
  set sample_id, reads from reformat_PE_in 
  
  output:
  set sample_id, "${sample_id}_phase1_input_PE.fasta.gz" into reformation_PE
  
  script:
""" 
reformat.sh in1=${reads[0]} in2=${reads[1]} out=${sample_id}_phase1_input_PE.fasta addslash
sleep 1s
gzip ${sample_id}_phase1_input_PE.fasta
sleep 1s
"""
}

process reformat_fastq_reads_S{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reformation", mode:'link'

  input:
  set sample_id, reads from reformat_S_in 
  
  output:
  set sample_id, "${sample_id}_phase1_input_S.fasta.gz" into reformation_S
  
  script:
""" 
reformat.sh in=${reads[2]} out=${sample_id}_phase1_input_S.fasta
sleep 1s
gzip ${sample_id}_phase1_input_S.fasta
sleep 1s
"""
}

process blastn_phase1_PE{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blast_results", mode:'link'

  input:
  set sample_id, reads from reformation_PE 
  
  output:
  set sample_id, "${sample_id}_blastn_phase1_PE.out" into blastn_phase1_PE
  
  script:
""" 
cp ${reads} ./tmp.fasta.gz
sleep 1s
gunzip tmp.fasta.gz
sleep 1s
blastn -db ${params.small_db_nucl} -query tmp.fasta -num_threads ${task.cpus} -outfmt 6 > "${sample_id}_blastn_phase1_PE.out"
sleep 1s
rm tmp.fasta
"""
}

process blastn_phase1_S{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blast_results", mode:'link'

  input:
  set sample_id, reads from reformation_S 
  
  output:
  set sample_id, "${sample_id}_blastn_phase1_S.out" into blastn_phase1_S
  
  script:
""" 
cp ${reads} ./tmp.fasta.gz
sleep 1s
gunzip tmp.fasta.gz
sleep 1s
blastn -db ${params.small_db_nucl} -query tmp.fasta -num_threads ${task.cpus} -outfmt 6 > "${sample_id}_blastn_phase1_S.out"
sleep 1s
rm tmp.fasta
"""
}

tmp = blastn_phase1_PE.combine(blastn_phase1_S, by: 0)
blastn_phase1_with_fq = tmp.combine(extraction_in, by: 0)

process extract_fastq_reads_for_phase1{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/phase1_fa_to_fq", mode:'link'

  input:
  set sample_id, blastn_PE, blastn_S, reads from blastn_phase1_with_fq 
  
  output:
  set sample_id, "${sample_id}_phase1_read_1.fq.gz", "${sample_id}_phase1_read_2.fq.gz", "${sample_id}_phase1_read_unpaired.fq.gz" into phase1_fq_reads
  script:
"""
cat ${blastn_PE} | awk '{ print \$1 }' | awk '!seen[\$0]++' > list_PE.txt
cat ${blastn_S} | awk '{ print \$1 }' | awk '!seen[\$0]++' > list_S.txt
seqtk subseq ${reads[0]} list_PE.txt | gzip > ${sample_id}_phase1_read_1.fq.gz
seqtk subseq ${reads[1]} list_PE.txt | gzip > ${sample_id}_phase1_read_2.fq.gz
seqtk subseq ${reads[2]} list_S.txt | gzip > ${sample_id}_phase1_read_unpaired.fq.gz
sleep 10
rm list_PE.txt list_S.txt
"""
}

phase1_fq_reads.into{asm_spades_in;
asm_megahit_in;
reads_for_mapping}

process asm_metaviralspades{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/assembly", mode:'link'

  input:
  set sample_id, read_1, read_2, read_un from asm_spades_in 
  
  output:
  set sample_id, "${sample_id}_asm_metaspades.fasta" into asm_spades_out
  file "spades_assembly" optional true into asm_spades_dir
  script:
""" 
\$assembler -1 ${read_1} -2 ${read_2} -s ${read_un} -t ${task.cpus} -m ${params.memory} -o spades_assembly
cp ./spades_assembly/contigs.fasta ./${sample_id}_asm_metaspades.fasta
"""
}

process asm_megahit{
  tag { "${sample_id}" }
  publishDir "${params.publish_base_dir}/${sample_id}/assembly", mode: 'link'

  input:
  set sample_id, read_1, read_2, read_un from asm_megahit_in

  output:
  set sample_id, "${sample_id}_asm_megahit.fasta" optional true into asm_megahit
  file "megahit_assembly" optional true into asm_megahit_dir

  script:
  """
  megahit -t ${task.cpus} --presets meta-sensitive -1 ${read_1} -2 ${read_2} -r ${read_un} --cpu-only  -o megahit_assembly
  if [ -s megahit_assembly/final.contigs.fa ]; then ln megahit_assembly/final.contigs.fa ${sample_id}_asm_megahit.fasta; fi
  """
}

process modify_asm_megahit{
  tag { "${sample_id}" }
  publishDir "${params.publish_base_dir}/${sample_id}/assembly", mode: 'link'

  input:
  set sample_id, asm from asm_megahit

  output:
  set sample_id, "${sample_id}_asm_megahit_modified.fasta" optional true into asm_megahit_mod_out

  script:
  """
  cat ${asm} | sed -e 's/\\s\\+/,/g' > "${sample_id}_asm_megahit_modified.fasta"
  """
}

asm_megahit_mod_out.into{asm_megahit_mod;
asm_megahit_comb}

asm_spades_out.into{asm_spades;
asm_spades_comb}

asm_for_mapping = asm_spades_comb.combine(asm_megahit_comb, by: 0)
map_reads_to_contigs = asm_for_mapping.combine(reads_for_mapping, by: 0)

process map_reads_to_contigs{
  tag { "${sample_id}" }

  publishDir "${params.publish_base_dir}/${sample_id}/assembly/reads_mapped_to_contigs", mode:'link'

  input:
  set sample_id, spades_contigs, megahit_contigs, read_PE_1, read_PE_2, read_single  from map_reads_to_contigs

  output:
  set sample_id, "${sample_id}_metaspades_reads_to_contigs.sam.gz", "${sample_id}_megahit_reads_to_contigs.sam.gz" into map_reads_to_contigs_out

  script:
  """
  bbwrap.sh ref=${spades_contigs} in=${read_PE_1},${read_single} in2=${read_PE_2},null out=${sample_id}_metaspades_reads_to_contigs.sam.gz usejni=t kfilter=22 subfilter=15 maxindel=80
  sleep 5s
  bbwrap.sh ref=${megahit_contigs} in=${read_PE_1},${read_single} in2=${read_PE_2},null out=${sample_id}_megahit_reads_to_contigs.sam.gz usejni=t kfilter=22 subfilter=15 maxindel=80
  """
}

map_reads_to_contigs_out.into{asm_mapping_stats_in;
asm_per_ctg_coverage_in;
parse_sam_mapping}

process asm_mapping_stats{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/assembly/reads_mapped_to_contigs", mode:'link'

  input:
  set sample_id, spades_mapping, megahit_mapping from asm_mapping_stats_in

  output:
  set sample_id, "${sample_id}_metaspades_flagstat.txt", "${sample_id}_megahit_flagstat.txt" into asm_mapping_stats_out

  script:
  """
  samtools flagstat ${spades_mapping} > ${sample_id}_metaspades_flagstat.txt
  sleep 5s
  samtools flagstat ${megahit_mapping} > ${sample_id}_megahit_flagstat.txt
  """
}

process asm_per_ctg_coverage{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/assembly/reads_mapped_to_contigs", mode:'link'

  input:
  set sample_id, spades_mapping, megahit_mapping from asm_per_ctg_coverage_in

  output:
  set sample_id, "${sample_id}_metaspades_reads_to_contigs.cov.txt", "${sample_id}_megahit_reads_to_contigs.cov.txt" into asm_per_ctg_coverage_out

  script:
  """
  pileup.sh in=${spades_mapping} out=${sample_id}_metaspades_reads_to_contigs.cov.txt 32bit=t
  sleep 5s
  pileup.sh in=${megahit_mapping} out=${sample_id}_megahit_reads_to_contigs.cov.txt 32bit=t
  """
}

process asm_bedtools_parse_reads_mapped_to_contigs{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/assembly/reads_mapped_to_contigs", mode:'link'

  input:
  set sample_id, spades_mapping, megahit_mapping from parse_sam_mapping

  output:
  set sample_id, "${sample_id}_bedtools_metaspades_reads_mapped_to_contigs.txt.gz", "${sample_id}_bedtools_megahit_reads_mapped_to_contigs.txt.gz" into parse_sam_mapping_out

  script:
  """
  bedtools bamtobed -i ${spades_mapping} | gzip > ${sample_id}_bedtools_metaspades_reads_mapped_to_contigs.txt.gz
  sleep 5s
  bedtools bamtobed -i ${megahit_mapping} | gzip > ${sample_id}_bedtools_megahit_reads_mapped_to_contigs.txt.gz
  """
}

process blastn_phase2_metaviralspades{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blast_results", mode:'link'

  input:
  set sample_id, reads from asm_spades 
  
  output:
  set sample_id, "${sample_id}_blastn_phase2_spades.out" into blastn_phase2_spades
  
  script:
""" 
blastn -db ${params.large_db_nucl} -query ${reads} -num_threads ${task.cpus} -outfmt 6 > "${sample_id}_blastn_phase2_spades.out"
"""
}

process blastn_phase2_megahit{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blast_results", mode:'link'

  input:
  set sample_id, reads from asm_megahit_mod 
  
  output:
  set sample_id, "${sample_id}_blastn_phase2_megahit.out" into blastn_phase2_megahit
  
  script:
""" 
blastn -db ${params.large_db_nucl} -query ${reads} -num_threads ${task.cpus} -outfmt 6 > "${sample_id}_blastn_phase2_megahit.out"
"""
}

blastn_phase2_results = blastn_phase2_spades.combine(blastn_phase2_megahit, by: 0)

process process_phase2_results{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blast_results", mode:'link'

  input:
  set sample_id, spades_res, megahit_res from blastn_phase2_results 
  
  output:
  set sample_id, "${sample_id}_blastn_phase2_spades_datamash_filt.txt", "${sample_id}_blastn_phase2_megahit_datamash_filt.txt", "${sample_id}_phase2_accession_num_list.txt" into phase2_processed_results
  
  script:
""" 
cat ${spades_res} | tr '.' ',' | datamash -g 1 min 11 -f | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$4 "\t" \$5 "\t" \$6 "\t" \$7 "\t" \$8 "\t" \$9 "\t" \$10 "\t" \$11 "\t" \$12}' | tr ',' '.' > ${sample_id}_blastn_phase2_spades_datamash_filt.txt

cat ${megahit_res} | tr '.' ',' | datamash -g 1 min 11 -f | awk '{print \$1 "\t" \$2 "\t" \$3 "\t" \$4 "\t" \$5 "\t" \$6 "\t" \$7 "\t" \$8 "\t" \$9 "\t" \$10 "\t" \$11 "\t" \$12}' | tr ',' '.' > ${sample_id}_blastn_phase2_megahit_datamash_filt.txt

cat ${sample_id}_blastn_phase2_megahit_datamash_filt.txt ${sample_id}_blastn_phase2_spades_datamash_filt.txt | awk '{print \$2}' | awk '!seen[\$0]++' > ${sample_id}_phase2_accession_num_list.txt

"""
}

phase2_processed_results.into{fetch_names_local_DB;
fetch_names_esearch;
initial_table}

process fetch_fasta_headers_local_DB{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/taxonomy/local_DB", mode:'link'

  input:
  set sample_id, spades_filt, megahit_filt, acc_list from fetch_names_local_DB
 
  output:
  set sample_id, "${sample_id}_organism_names_local_DB.tsv", "${sample_id}_acc_nr_not_found.txt", "${sample_id}_acc_nr_exist_for_esearch.txt" into local_DB_organism_name_fetched
  script:
"""
count=1
echo "Accession_nr	Title	Organism" > tmp.tsv
echo "Accession_numbers_not_found_in_nr_db" > acc_nr_not_found.txt
echo "Existing acession number for esearch" > acc_nr_exist_for_esearch.txt
for line in \$(cat ${acc_list})
do
	title_blast=\$(blastdbcmd -entry \$line -db ${params.large_db_nucl} -range 1-1)
	title=\$(echo \$title_blast|awk -F \$line '{print \$2}'| awk -F '>' '{print \$1}')
	if [ -z "\$title" ]
	then
		if [ -z "\$title_blast" ]
		then
			echo "\$line" >> acc_nr_not_found.txt
		else
			echo "\$line" >> acc_nr_exist_for_esearch.txt
		fi
	else
		organism=\$(echo "\$title"|awk -F '[' '{print \$2}'|awk -F ']' '{print \$1}')
		organism=\${organism// /_}
		title=\${title// /_}
		title=\${title//:1-1_/}
		echo \$line	\$title	\$organism >>tmp.tsv
	fi
done
tr ' ' '\t' <tmp.tsv >"${sample_id}_organism_names_local_DB.tsv"
rm tmp.tsv
mv acc_nr_not_found.txt ${sample_id}_acc_nr_not_found.txt
mv acc_nr_exist_for_esearch.txt ${sample_id}_acc_nr_exist_for_esearch.txt

"""
}

process fetch_taxonomy_esearch{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/taxonomy/esearch", mode:'link'

  input:
  set sample_id, spades_filt, megahit_filt, acc_list from fetch_names_esearch
 
  output:
  set sample_id, "${sample_id}_phase2_esearch_taxonomy.tsv" into esearch_taxonomy_fetched
  script:
"""
API_KEY="${params.api_key}"
touch test.tsv
echo "Accession_nr	Title	Organism	Division	Rank" > header.txt
for line in \$(cat ${acc_list})
do
	esearch  -db nucleotide -query \$line -api_key=\$API_KEY|esummary > temp.txt
	sleep 1s
	title=\$(cat temp.txt | grep Title| awk -F '</Title>' '{print \$1}' | awk -F '<Title>' '{print \$2}')
	organism=\$(cat temp.txt | grep Organism| awk -F '</Organism>' '{print \$1}'|awk -F '<Organism>' '{print \$2}')
	esearch  -db taxonomy -query "\$organism" -api_key=\$API_KEY|esummary > temp.txt
	sleep 1s
	rank=\$(cat temp.txt | grep Rank| awk -F '</Rank>' '{print \$1}' | awk -F '<Rank>' '{print \$2}')
	div=\$(cat temp.txt | grep Division| awk -F '</Division>' '{print \$1}'|awk -F '<Division>' '{print \$2}')
	organism=\${organism// /_}
	title=\${title// /_}
	rank=\${rank// /_}
	div=\${div// /_}
	echo -e \$line"\t"\$title"\t"\$organism"\t"\$div"\t"\$rank >> test.tsv
	sleep 3s   
done
cat header.txt test.tsv > "${sample_id}_phase2_esearch_taxonomy.tsv"
rm temp.txt test.tsv header.txt

"""
}

classification_data = esearch_taxonomy_fetched.combine(initial_table, by: 0)

process preprocess_classification_data{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/taxonomy/classification", mode:'link'

  input:
  set sample_id, taxonomy, spades, megahit, acc_list from classification_data
 
  output:
  set sample_id, "${sample_id}_both_assemblers_blastn_results_sorted_by_acc_num.txt", "${sample_id}_phase2_esearch_taxonomy_sorted.tsv" into preprocessed_classification
  script:
"""
cat ${spades} ${megahit} > tmp.txt
cat tmp.txt | sort -k2 > "${sample_id}_both_assemblers_blastn_results_sorted_by_acc_num.txt"
cat ${taxonomy} | sort -k1 > "${sample_id}_phase2_esearch_taxonomy_sorted.tsv"
rm tmp.txt
"""
}

process classified_asm_results{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/tables/tsv_tables", mode:'link'

  input:
  set sample_id, asm_res, tax_list from preprocessed_classification
 
  output:
  set sample_id, "${sample_id}_all_asm_results.tsv", "${sample_id}_viral_asm_results.tsv" into classified_asm
  script:
"""
#!/usr/bin/python3
f1 = open("${tax_list}","r+")
taxonomy = f1.readlines()
f2 = open("${asm_res}","r+")
asm_results = f2.readlines()
f1.close()
f2.close()
header = "sample_id\\tseq_id\\tref_id\\ttitle\\tscientific_name\\trank\\tdivison\\tperc_of_identical_matches\\tevalue\\tmismatches\\talignment_len\\tseq_len\\tassembler\\n"
f1 = open("${sample_id}_all_asm_results.tsv","a+")
f2 = open("${sample_id}_viral_asm_results.tsv","a+")
f1.write(header)
f2.write(header)
c = ['-'] * 13 #container
for tax in taxonomy:
    for asm_r in asm_results:
        tmp = asm_r.split()
        tmp2 = tax.split()
        if str(tmp2[0]) == str(tmp[1]):
            c[0] = "${sample_id}"
            if "node" in tmp[0].lower(): #adding sequence id to first column
                c[1] = tmp[0].split("_len")[0]
                # sequence len
                c[11] = tmp[0].split("_length_")[1].split("_cov")[0]
                # assembler
                c[12] = "metaviralspades"
            else:
                c[1] = tmp[0].split(".flag")[0]
                # sequence len
                c[11] = tmp[0].split(".len=")[1]
                # assembler
                c[12] = "megahit"
           # adding reference id to 2nd column
            c[2] = tmp[1]
            # title
            c[3] = tmp2[1]
            # scientific name
            c[4] = tmp2[2]      
            if len(tmp2) == 5:
                # rank
                c[5] = tmp2[4]
                # divison
                c[6] = tmp2[3] 
            # identical
            c[7] = tmp[2]
            c[8] = tmp[10]
            # mismatches
            c[9] = tmp[4] 
            # alignment length
            c[10] = tmp[3]
            str_tax = str(c[0])+"\\t"+str(c[1])+"\\t"+str(c[2])+"\\t"+str(c[3])+"\\t"+str(c[4])+"\\t"+str(c[5])+"\\t"+str(c[6])+"\\t"+str(c[7])+"\\t"+str(c[8])+"\\t"+str(c[9])+"\\t"+str(c[10])+"\\t"+str(c[11])+"\\t"+str(c[12])+"\\n"
            f1.write(str_tax)
            str_cont=""
            str_cont = str_cont.join(c[:-1])
            if "virus" in str_cont.lower() or "viral" in str_cont.lower() or "phage" in str_cont.lower():
                f2.write(str_tax)
            c = ['-'] * 13
f1.close()
f2.close()
"""
}

classified_asm.into{tsv_into_html_each;
combine_viral_tsv;
combine_all_tsv}

html_data_each = tsv_into_html_each.combine(html_each)

process TSV_to_HTML_each_sample{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/tables/html_tables", mode:'link'

  input:
  set sample_id, all, viral, hvalue, html from html_data_each
 
  output:
  set sample_id, "${sample_id}_all_asm_results.html", "${sample_id}_viral_asm_results.html" into html_each_out
  script:
"""
#!/usr/bin/python3
import string
import csv
import json
import collections

OrderedDict = collections.OrderedDict

def HTML_table(src, header, htmlEndStr, name):
    
    f = open(src)
    lines = f.readlines()
    f.close()
    f1 = open(header)
    header_lines = f1.readlines()
    f2 = open(htmlEndStr)
    end_lines = f2.readlines()
    f1.close()
    f2.close()
    
    with open(name, 'w') as f:
        for h_line in header_lines:
            f.write(h_line)
        for line in lines:
            f.write(line)
        for e_line in end_lines:
            f.write(e_line)
    f.close()
    return None

def TSV_file_into_JSON(src, dst, header):
    data = []
    with open(src, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\\t', quotechar='"')
        for row in reader:
            if row[0] == 'sample_id':
                print("\\n")
            else:
                if row[0].strip()[0] == '#':  #
                    continue
                row = filter(None, row)
                data.append(OrderedDict(zip(header, row)))

    with open(dst, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=2)
    return None

header = ['sample_id', 'seq_id', 'ref_id', 'title', 'scientific_name', 'rank', 'divison', 'perc_of_identical_matches', 'evalue', 'mismatches', 'alignment_len', 'seq_len', 'assembler']
src = "${all}"
dst = "all.json"  
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${sample_id}_all_asm_results.html")

src = "${viral}"
dst = "viral.json"  
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${sample_id}_viral_asm_results.html")

"""
}

process collect_all_tsv_tables{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tables/tsv_tables", mode:'link'

  input:
  file "tsv_table" from combine_all_tsv.map{it[1]}.collect()
 
  output:
  file "${params.project_id}_all_asm_results.tsv" into all_tsv_collected
  script:
"""
touch tmp.tsv
echo -e "sample_id\\tseq_id\\tref_id\\ttitle\\tscientific_name\\trank\\tdivison\\tperc_of_identical_matches\\tevalue\\tmismatches\\talignment_len\\tseq_len\\tassembler" > header.tsv

for sample_file in ${tsv_table}
do
	cat \$sample_file | sed '1d' >> tmp.tsv
done

cat tmp.tsv | sort -k1,1 -k3,3 > tmp2.tsv
cat header.tsv tmp2.tsv > "${params.project_id}_all_asm_results.tsv"
rm tmp.tsv tmp2.tsv header.tsv
"""
}

process collect_viral_tsv_tables{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tables/tsv_tables", mode:'link'

  input:
  file "tsv_table" from combine_viral_tsv.map{it[2]}.collect()
 
  output:
  file "${params.project_id}_viral_asm_results.tsv" into viral_tsv_collected
  script:
"""
#!/bin/bash
touch tmp.tsv
echo -e "sample_id\\tseq_id\\tref_id\\ttitle\\tscientific_name\\trank\\tdivison\\tperc_of_identical_matches\\tevalue\\tmismatches\\talignment_len\\tseq_len\\tassembler" > header.tsv

for sample_file in ${tsv_table}
do
	cat \$sample_file | sed '1d' >> tmp.tsv
done

cat tmp.tsv | sort -k1,1 -k3,3 > tmp2.tsv
cat header.tsv tmp2.tsv > "${params.project_id}_viral_asm_results.tsv"
rm tmp.tsv tmp2.tsv header.tsv
"""
}

tsv_collection_into_html = all_tsv_collected.combine(viral_tsv_collected).combine(html_all)

process TSV_to_HTML_all_collected{
  tag {"All"}

  publishDir "${params.publish_base_dir}/All/tables/html_tables", mode:'link'

  input:
  set "all", "viral", "hvalue", "html" from tsv_collection_into_html
 
  output:
  set "${params.project_id}_all_asm_results.html", "${params.project_id}_viral_asm_results.html" into html_all_out
  script:
"""
#!/usr/bin/python3
import string
import csv
import json
import collections

OrderedDict = collections.OrderedDict

def HTML_table(src, header, htmlEndStr, name):
    
    f = open(src)
    lines = f.readlines()
    f.close()
    f1 = open(header)
    header_lines = f1.readlines()
    f2 = open(htmlEndStr)
    end_lines = f2.readlines()
    f1.close()
    f2.close()
    
    with open(name, 'w') as f:
        for h_line in header_lines:
            f.write(h_line)
        for line in lines:
            f.write(line)
        for e_line in end_lines:
            f.write(e_line)
    f.close()
    return None

def TSV_file_into_JSON(src, dst, header):
    data = []
    with open(src, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\\t', quotechar='"')
        for row in reader:
            if row[0] == 'sample_id':
                print("\\n")
            else:
                if row[0].strip()[0] == '#':  #
                    continue
                row = filter(None, row)
                data.append(OrderedDict(zip(header, row)))

    with open(dst, 'w') as jsonfile:
        json.dump(data, jsonfile, indent=2)
    return None

header = ['sample_id', 'seq_id', 'ref_id', 'title', 'scientific_name', 'rank', 'divison', 'perc_of_identical_matches', 'evalue', 'mismatches', 'alignment_len', 'seq_len', 'assembler']
src = "${all}"
dst = "all.json"  
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_all_asm_results.html")

src = "${viral}"
dst = "viral.json"  
TSV_file_into_JSON(src, dst, header)
HTML_table(dst, "${html[1]}", "${html[0]}", "${params.project_id}_viral_asm_results.html")

"""
}



