# variantSeeker

A pipeline that is designed to find variants / types based on two phases using reference databases to target the most common variants of interest. Phase 1 collects all reads that match with the phase 1 database. The database for phase 1 contains only the sequences we are intrested in and is a much smaller database than phase 2. The reads collect in phase 1 will be used in two metagenomic assemblers and the contigs from the metagenomic assembly will go into phase 2 where we are using 'nt' (Partially non-redundant nucleotide sequences) database. Based on the results from phase 2 the contigs will be classified and presented in a html and tsv tables.

![alt text](/IMG/UML_diagrams/variantSeeker.png)


The UML (Unified Modeling Language) diagram above displays the steps in the pipeline. "Input data" is the input reads and the data should be trimmed, quality checked and human DNA/reads should be filtered out before running the pipeline. Here we can see a summary of all the steps.

Below we can see an example of a table in the last process and what we were interested in were TTV (torque teno virus) viruses.
![alt text](/IMG/tables/variantSeeker_html_table.png)

## Software requirements 
 Software required to run all processes in the pipeline.
 - [Nextflow DSL1](https://www.nextflow.io/)
 - [Python3](https://www.python.org/downloads/)
 - [megahit](https://github.com/voutcn/megahit)
 - [metaviralspades](https://cab.spbu.ru/software/spades/)
 - [seqtk](https://github.com/lh3/seqtk)
 - [BBmap tools](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/)
    - reformat.sh
    - bbwrap.sh
    - pileup.sh
 - [samtools](http://www.htslib.org/)
 - [datamash](https://www.gnu.org/software/datamash/)
 - [Entrez Direct: E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
    - esearch
 - [BLAST blastdbcmd](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
 - [BLAST blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

#### Software versions currently tested on
| Software   | Version |
| --------   | ------- |
| BBmap      | 38.68   |
| BLAST      | 2.9.0+  |
| Seqtk      | 1.3-r106|
| metaSPAdes | v3.14.1 |
| MEGAHIT    | v1.2.9  |
| samtools   | 1.9     |
| bedtools   | v2.29.0 |
| GNU datamash | 1.1.0 |
| entrez-direct | 13.9 |
| Python3    | 3.7.8   |

|Python packages|
| ------------- |
| string        |
| csv           |
| json          |
| collections   |

## Database requirements 

You will need to download several databses to be able to run the pipeline. It may be a good idea to create a folder for each database.

#### 1. Phase 1 databse 
For phase 1, we will build a database where we use sequences that belong to a variant / type that we are interested in finding among our input samples. In this example below, we are looking for TTV viruses. You will need to replace "<reference.fasta>" with the fasta file which contains your sequences and name the database. Here I have named it "TTV_nucleotide".

```
makeblastdb -in <reference.fasta> -dbtype nucl -parse_seqids  -out TTV_nucleotide -title "TTV_nucleotide"
```

#### 2. Phase 2 database
For phase 2 we will need the entire "nt" (Partially non-redundant nucleotide sequences) database. First we need to download fasta file and then run "makeblastdb" as we did earlier in phase 1.

Link to "nt" database:
```
https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz
```
How to create databse for phase 2:

```
makeblastdb -in nt.fasta -dbtype nucl -parse_seqids  -out nt -title "nt"
```

#### 3. blastdbcmd
Here we only need to reuse the database "nt" we obtained in phase 2.

## NCBI API KEY
It's easy to obtain a key. All you have to do is create an account at NCBI, log in and generate a key.

First, go to the website below.
```
https://www.ncbi.nlm.nih.gov/
```
Then click on "Sign in to NCBI" at the top right corner. Once you have logged in / registered, click on your username at the top right corner. You will be forwarded to your account where you will scroll down to "API Key Management" and click on create an API Key. You will then need to paste the key into the configuration file "variantSeeker.config" which you can read more about in "Configuration parameters".

## Configuration file
In our configuration file "variantSeeker.config" we can add a profile. Each profile can be customized according to the available resources of your system. You can create your own profile through the following example below which you can find in "variantSeeker.config" file.

```
profiles {

  amanj {
    includeConfig 'conf/amanj.variantSeeker.config'
  }
  
  othello {
    includeConfig 'conf/othello.variantSeeker.config'
  }
}

```
Above we can see two profiles. One of them is included in the repository created for "amanj". You can follow the exemples below on how to modify "amanj" and create youre own file and include this later in "variantSeeker.config" file.

### Configuration parameters
In the configuration file "variantSeeker.config" we also have parameters that are used during the pipeline run.

```
/* Pipeline running parameters */
params{
  small_db_nucl='TTV_nucleotide'
  large_db_nucl='nt'
  publish_base_dir='variantSeeker_output'
  fastq_dir='input_reads'
  html_dir='input_html'
  memory=110
  project_id='P13408'
  api_key="API KEY"
}
```
"small_db_nucl" and "large_db_nucl" should be the names you chose for the databases in phase 1 and phase 2 (see database requirements). Remember that phase 1 is smaller than phase 2.

"project_id" is the label for your project.

"api_key" is required for fetching taxonomic data from an external database connected to the E-utilities and you should put the key in "api_key" parameter. The API KEY allows us to increase the number of request. See "NCBI API KEY" on how to obtain a key. 

### Database pathways 
In our configuration file found in conf/amanj.variantSeeker.config we will need to add all the paths for the databases.
Below you can see files and folders you will need to add full path to.

For "blastn_phase1_PE" and "blastn_phase1_S" you will need the complete path to phase 1 database folder and the other three requires a complete path to the folder where the phase 2 database is located. Note that you will need to add the database names to the configuration parameters which you can read more about in "Configuration parameters".

```
    withName: blastn_phase1_PE{
       beforeScript='export BLASTDB=/FullPathToFolder/database/small_DB_Anelloviruses_nucleotide'
       cpus=25
    }
    
    withName: blastn_phase1_S{
       beforeScript='export BLASTDB=/FullPathToFolder/database/small_DB_Anelloviruses_nucleotide'
       cpus=14
    }
    withName: blastn_phase2_megahit{
       beforeScript='export BLASTDB=/FullPathToFolder/Database/nt'
       cpus=14
    }

    withName: blastn_phase2_metaviralspades{
       beforeScript='export BLASTDB=/FullPathToFolder/Database/nt'
       cpus=14
    }

    withName: fetch_fasta_headers_local_DB{
       beforeScript='export BLASTDB=/FullPathToFolder/Database/nt'
    }

```

### Software pathways 
In our configuration file found in conf/amanj.discovery.config we will need to add all the paths for all softwares if they are not included in your PATH variable or conda environment.

Below you can see an example of a process which requires a path the BBmap tools. The default value for number of threads are 1 in each process if you want to increase you will need to add "cpus = 8" to increase it and in this case it is increased to 8 threads. You will need to do this for all processes if they are not part of your current PATH variable or conda environment. All processes are executed as a child process and adding the line "beforeScript='export PATH="/PathToFolder/tools/BBMap/38.68/bbmap:$PATH"'" does not save it to your PATH variable. It is only temprorary during the pipeline run. 
```
    withName: asm_map_reads_to_contigs{
		beforeScript='export PATH="/PathToFolder/tools/BBMap/38.68/bbmap:$PATH"'
        cpus = 8
    }
```
#### metaspades pathway 
For the metaspades assembler I made a variable that requries a complete path to the "metaspades.py" program in "conf/amanj.variantSeeker.config" as exemplified below.
```
    withName: asm_metaviralspades{
       beforeScript='export assembler=/FullPathToFile/metaviralspades/spades/bin/metaspades.py'
       cpus=10
    }
```
## Running pipeline
The user should create a folder one called 'input_reads' and store all samples there. Each sample should have a folder labeled with sample name / samle ID. Each sample folder should contain paired-end fastq compressed gzip (GNU Zip) files and a file with unpaired reads in fastq gzip format. Each fastq file need to be labled with same sample names but different file extension names "_1.fq.gz", "_2.fq.gz" and "_unpaired.fq.gz" for the pipeline to recognize the sample.

If you prefer "fastq.gz" over "fq.gz" you can change in the line below and this line can be found in begning of "variantSeeker.nf" file.

```
fastq_files = Channel.fromFilePairs("${params.fastq_dir}/**/*_{1,2,unpaired}.fq.gz",size:3)
```

To run the pipeline in command line:
```
nextflow -C variantSeeker.config run variantSeeker.nr -profile amanj
```
To run the pipeline in command line and resume from cache memory:
```
nextflow -C variantSeeker.config run variantSeeker.nr -profile amanj -resume
```





