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
    #1 bootstrap for now to test
        "kallisto quant -i {input.index} -o kallisto/{wildcards.sample} -b 1 -t 1 {input.r1} {input.r2}"
        
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

# rule copy_hello_world:
#     input:
#         input_file="/home/slarosa/snakemake_demo/input.txt"
#     output:
#         output_file="/home/slarosa/snakemake_demo/output_dir/hello_world.txt"
#     shell:
#         "cp {input.input_file} {output.output_file}"

# rule clean:
#     shell:
#         "rm -r ./output_dir"


rule bowtie_build:
    input:
        fasta="GCA_000845245.1_ViralProj14559_genomic.fna"
    shell:
        "bowtie2-build {fasta} HCMV"

rule bowtie_map:
    input:
        fq1="fastq_files/{sample}_R1.fastq",
        fq2="fastq_files/{sample}_R2.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bowtie2 --quiet -x HCMV -1 {input.fq1} -2 {input.fq2} -S {output} --al-conc-gz {sample}_mapped_%.fq.gz"
