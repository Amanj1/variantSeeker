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
