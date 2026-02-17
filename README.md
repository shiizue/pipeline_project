# pipeline_project
COMP438 independent project
### 1
created txt file with the accession numbers, one per line
used a script to fasterq-dump each one in a loop (Googled 'fasterq dump multiple accessions)

### 2
created a python script (cds.py) to use biopython to get genome record (biopython docs)
get refseq protein id and cds from record and write to a new fasta and pipeline report
(https://github.com/peterjc/biopython_workshop/blob/master/using_seqfeatures/README.rst)
wrote rules in snakefile to create index with kallisto

### 3 
wrote rule in snakefile to quantify the TPM of each CDS from the fasta file from step 2
create r script to use sleuth to compare the expressed genes between the 2 conditions, which outputs a tab-delimited table (in-class code)
added rule in snakefile to call the R script

converted script from step 1 into rule in snakemake
created rule in snakemake to write all the results from each step into PipelineReport