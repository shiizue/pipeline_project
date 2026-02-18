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
os.makedirs("assembly", exist_ok=True)

#run all rules until we get the final PipelineReport.txt
rule all:
    input:
        "PipelineReport.txt"

#rule to write PipelineReport by calling write_report.py
rule write_report:
    input:
    #include other steps later
        cds="results/cds.txt",
        sleuth="results/sleuth_results.txt",
        bowtie="results/bowtie_log.txt",
        #doesn't actually write anything from SPAdes, just needs to know it has to run this rule
        assembly=expand("assembly/{sample}/contigs.fasta", sample=samples),
        blast="results/blast_results.txt"
    output:
        "PipelineReport.txt"
    run:
        with open(output[0], "w") as f:
            for input_file in [input.cds,input.sleuth,input.bowtie,input.blast]:
                with open(input_file) as r:
                    f.write(r.read())
                f.write("\n")


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
#this is what the next step is going to map to
rule bowtie_build:
    input:
        fasta="data/GCA_000845245.1_ViralProj14559_genomic.fna"
    output:
        expand("bowtie/HCMV.{ext}", ext=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
    shell:
        "bowtie2-build {input.fasta} bowtie/HCMV"

#keep only reads that map to the HCMV genome
rule bowtie_map:
    input:
        fq1="fastq_files/{sample}_1.fastq",
        fq2="fastq_files/{sample}_2.fastq",
        index=expand("bowtie/HCMV.{ext}", ext=["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"])
    output:
        mapped_r1="mapped_reads/{sample}_mapped_1.fastq",
        mapped_r2="mapped_reads/{sample}_mapped_2.fastq",
        mapped_log="mapped_reads/{sample}_bowtie_log.txt"
    shell:
    #just want filtered reads in a fastq file, don't care for sam output
        "bowtie2 --quiet -x bowtie/HCMV -1 {input.fq1} -2 {input.fq2} --al-conc mapped_reads/{wildcards.sample}_mapped_%.fastq -S /dev/null 2> {output.mapped_log}"

rule bowtie_log:
    input:
        original=expand("fastq_files/{sample}_1.fastq", sample=samples),
        logs=expand("mapped_reads/{sample}_bowtie_log.txt", sample=samples)
    output:
        "results/bowtie_log.txt"
    run:
        with open(output[0], "w") as f:
            for sample in samples:
                # count before mapping
                with open(f"fastq_files/{sample}_1.fastq") as fastq:
                #reads only on 4th line
                    original_count = sum(1 for line in fastq) // 4

                # parse concordantly aligned reads from bowtie log
                #these are the reads that actually mapped and we're filtering for
                mapped = 0
                with open(f"mapped_reads/{sample}_bowtie_log.txt") as log:
                    for line in log:
                        if "concordantly exactly 1 time" in line or "concordantly >1 times" in line:
                            mapped += int(line.strip().split()[0])

                f.write(
                    f"Sample {sample} had {original_count} read pairs before and {mapped} read pairs after Bowtie2 filtering.\n")

rule spades:
    input:
        r1="mapped_reads/{sample}_mapped_1.fastq",
        r2="mapped_reads/{sample}_mapped_2.fastq"
    output:
        "assembly/{sample}/contigs.fasta"
    shell:
        "spades.py -k 127 -t 1 --only-assembler -1 {input.r1} -2 {input.r2} -o assembly/{wildcards.sample}/"

rule blast:
    input:
        "assembly/{sample}/contigs.fasta"
    output:
        "results/{sample}_blast.txt"
    script:
        "scripts/blast.py"

rule blast_results:
    input:
        expand("results/{sample}_blast.txt",sample=samples)
    output:
        "results/blast_results.txt"
    run:
        with open(output[0], "w") as f:
            for sample, blast_file in zip(samples, input):
                f.write(f"{sample}:\n")
                f.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
                with open(blast_file) as b:
                    f.write(b.read())
                f.write("\n")