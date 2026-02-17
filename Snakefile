samples=['SRR5660030','SRR5660033','SRR5660044','SRR5660045']

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
        "kallisto quant -i {input.index} -o {output} -b 10 -t 2 {input.r1} {input.r2}"
        
#rule to run the R script to use sleuth to compare the 2 conditions
rule sleuth:
    input:
    #take every abundance.h5 file for each sample
        expand("kallisto/{sample}/abundance.h5", sample=samples),
        "sleuth_table.txt"
    output:
    #this result will get written to the final PipelineReport.txt
        "results/sleuth_results.txt"
    shell:
        "Rscript scripts/sleuth.R"
