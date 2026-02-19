import os
from Bio import SeqIO

contigs = snakemake.input.contigs
sample = snakemake.wildcards.sample
outfile = snakemake.output[0]
path = f"blast/betaherpesvirinae"

records = list(SeqIO.parse(contigs, "fasta"))
longest_contig = max(records, key=lambda r: len(r.seq))

longest_contig_fasta = f"assembly/{sample}/longest_contig.fasta"
with open(longest_contig_fasta, "w") as f:
    SeqIO.write(longest_contig, f, "fasta")

run_blastn = f'blastn -query {longest_contig_fasta} -db {path} -out {outfile} -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_hsps 1 -num_alignments 5'

os.system(run_blastn)
