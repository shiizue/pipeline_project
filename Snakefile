samples=['SRR5660030','SRR5660033','SRR5660044','SRR5660045']

rule all:
    input:
        "PipelineReport.txt"

rule write_report:
    input:
        cds="results/cds.txt",
        sleuth="results/sleuth_results.txt"
    output:
        "PipelineReport.txt"
    script:
        "scripts/write_report.py"

rule fasterq_dump:
    output:
        "fastq_files/{sample}_1.fastq",
        "fastq_files/{sample}_2.fastq"
    shell:
        "fasterq-dump {wildcards.sample} --outdir fastq_files/"

rule get_cds:
    output:
        fasta = "results/hcmv_cds.fasta",
        cds="results/cds.txt"
    script:
        "scripts/cds.py"  

rule kallisto_index:
    input:
        "results/hcmv_cds.fasta"
    output:
        "kallisto/hcmv.idx"
    shell:
        "kallisto index -i {output} {input}"

rule kallisto_quant:
    input:
        index = "kallisto/hcmv.idx",
        r1 = "fastq_files/{sample}_1.fastq",
        r2 = "fastq_files/{sample}_2.fastq"
    output:
        "kallisto/{sample}/abundance.h5"
    shell:
        "kallisto quant -i {input.index} -o {output} -b 10 -t 2 {input.r1} {input.r2}"
        
rule sleuth:
    input:
        expand("kallisto/{sample}/abundance.h5", sample=samples),
        "sleuth_table.txt"
    output:
        "results/sleuth_results.txt"
    shell:
        "Rscript scripts/sleuth.R"
