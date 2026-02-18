import os

samples=['SRR5660030','SRR5660033','SRR5660044','SRR5660045']
#dictionary to create table for sleuth later
conditions={'SRR5660030': '2dpi',
    'SRR5660033': '6dpi',
    'SRR5660044': '2dpi',
    'SRR5660045': '6dpi'
}

#set up folders for outputs
os.makedirs("results", exist_ok=True)
os.makedirs("mapped_reads", exist_ok=True)
os.makedirs("kallisto", exist_ok=True)
os.makedirs("data", exist_ok=True)

#run all rules until we get the final PipelineReport.txt
rule all:
    input:
        "PipelineReport.txt"

#rule to write PipelineReport by calling write_report.py
rule write_report:
    input:
    #include other steps later
        cds="results/cds.txt",
        sleuth="results/sleuth_results.txt"
    output:
        "PipelineReport.txt"
    script:
        "scripts/write_report.py"

#rule to run fasterq-dump for each fastq file by calling the above accession numbers
#puts fastq files in their own folder
rule fasterq_dump:
    output:
        "fastq_files/{sample}_1.fastq",
        "fastq_files/{sample}_2.fastq"
    shell:
        "fasterq-dump {wildcards.sample} --outdir fastq_files/"

#rule to get the total # of CDS and get a fasta file with protein_id and CDS's by calling cds.py
#this will get written to the final PipelineReport.txt
rule get_cds:
    output:
        fasta = "results/hcmv_cds.fasta",
        cds="results/cds.txt"
    script:
        "scripts/cds.py"  

#rule to use kallisto to build hcmv index using the fasta of the CDS's 
rule kallisto_index:
    input:
        "results/hcmv_cds.fasta"
    output:
        "kallisto/hcmv.idx"
    shell:
        "kallisto index -i {output} {input}"

#rule to quantify TPM using kallisto, referencing code from in class
rule kallisto_quant:
    input:
        index = "kallisto/hcmv.idx",
        r1 = "fastq_files/{sample}_1.fastq",
        r2 = "fastq_files/{sample}_2.fastq"
    output:
    #saves each output to a folder for each sample
        "kallisto/{sample}/abundance.h5"
    shell:
    #1 bootstrap for now to test
        "kallisto quant -i {input.index} -o kallisto/{wildcards.sample} -b 1 -t 1 {input.r1} {input.r2}"
        

#rule to make table for sleuth to read
rule sleuth_input_table:
    input:
        expand("kallisto/{sample}/abundance.h5", sample=samples)
    output:
        "data/sleuth_table.txt"
    run:
    #loop to write a tab delimited header and columns
        with open(output[0],"w") as f:
            f.write("sample\tpath\tcondition\n")
            for sample in samples:
                path =f'kallisto/{sample}'
                condition=conditions[sample]
                f.write(f'{sample}\t{path}\t{condition}\n')

#rule to run the R script to use sleuth to compare the 2 conditions
rule sleuth:
    input:
    #take every abundance.h5 file for each sample
        expand("kallisto/{sample}/abundance.h5", sample=samples),
        table="data/sleuth_table.txt"
    output:
    #this result will get written to the final PipelineReport.txt
        "results/sleuth_results.txt"
    shell:
        "Rscript scripts/sleuth.R {input.table}"

#build bowtie2 index from the HCMV genome fasta
rule bowtie_build:
    input:
        fasta="data/GCA_000845245.1_ViralProj14559_genomic.fna"
    output:
        expand("bowtie/HCMV.{ext}", ext=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
    shell:
        "bowtie2-build {input.fasta} bowtie/HCMV"

#keep only reads that map
rule bowtie_map:
    input:
        fq1="fastq_files/{sample}_1.fastq",
        fq2="fastq_files/{sample}_2.fastq",
        index=expand("bowtie/HCMV.{ext}", ext=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
    output:
        mapped_r1="mapped_reads/{sample}_mapped_1.fastq",
        mapped_r2="mapped_reads/{sample}_mapped_2.fastq",
    shell:
        "bowtie2 --quiet -x bowtie/HCMV -1 {input.fq1} -2 {input.fq2} -S mapped_reads/{wildcards.sample}_bowtie.sam --al-conc mapped_reads/{wildcards.sample}_mapped_%.fastq"
