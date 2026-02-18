import os
from Bio import SeqIO

download_data='datasets download virus genome taxon Betaherpesvirinae --refseq --include genome'
make_blast_db='makeblastdb -in betaherpesvirinae.fasta -out blast/betaherpesvirinae -title betaherpesvirina -dbtype nucl'

os.system(download_data)
os.system(make_blast_db)

contigs=snakemake.input.contigs
sample=snakemake.wildcards.sample
outfile=snakemake.output.blast

records = list(SeqIO.parse(contigs, "fasta"))
longest_contig=longest = max(records, key=lambda r: len(r.seq))

longest_contig_fasta=f"assembly/{sample}/longest_contig.fasta"
with open(longest_contig_fasta,"w") as f:
    SeqIO.write(longest_contig,f,"fasta")

run_blastn=f'blastn -query {longest_contig_fasta} -db betaherpesvirinae -out {outfile} --outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_hsps 1 -num_alignments 5'

os.system(run_blastn)