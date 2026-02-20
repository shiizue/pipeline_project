import os
from Bio import SeqIO

contigs = snakemake.input.contigs
sample = snakemake.wildcards.sample
outfile = snakemake.output[0]
path = f"blast/betaherpesvirinae"

#sort contigs by length to retrieve only the longest contig
records = list(SeqIO.parse(contigs, "fasta"))
longest_contig = max(records, key=lambda r: len(r.seq))

#write that contig to its own fasta file so we can query it into blast
longest_contig_fasta = f"assembly/{sample}/longest_contig.fasta"
with open(longest_contig_fasta, "w") as f:
    SeqIO.write(longest_contig, f, "fasta")

# get Subject accession, Percent identity, Alignment length, Start of alignment in query, End of alignment in query, Start of alignment in subject, End of alignment in subject, Bit score, E-value, and Subject Title
# only keep best alignment (hsp)
# only get top 5 hits
run_blastn = f'blastn -query {longest_contig_fasta} -db {path} -out {outfile} -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_hsps 1 -num_alignments 5'

os.system(run_blastn)
