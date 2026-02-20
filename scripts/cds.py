"""
cds.py
Uses BioPython to retrieve coding sequence (CDS) features from the HCMV genome.
Writes the CDS and their RefSeq protein IDs to a FASTA file for use in kallisto.
Counts the total CDS to be written to the final pipeline report.
"""

from Bio import Entrez, SeqIO

# if running locally, replace with your own email (but not really necessary)
Entrez.email = "skyel1005@gmail.com"
# I got the RefSeq from the NCBI page for GCF_000845245.1/ViralProj14559
refseq = "NC_006273.2"

# retrieve refseq data from NCBI
handle = Entrez.efetch(db="nucleotide", id=refseq, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")

cds_list = []
total_cds = 0

# loop through the "features", which includes the cds
# for each cds, add it to the total count, and also get its protein_id and its sequence
for feature in record.features:
    if feature.type == "CDS":
        total_cds += 1
        protein_id = feature.qualifiers["protein_id"][0]
        cds_seq = feature.extract(record.seq)
        cds_list.append([protein_id, cds_seq])

# create a fasta file with the previously retrieved protein_ids as headers for each cds and their sequences as the content
with open("results/hcmv_cds.fasta", "w") as f:
    for protein_id, cds_seq in cds_list:
        f.write(f">{protein_id}\n{cds_seq}\n")

# write number of cds to a txt file which will be written to the final PipelineReport.txt at the end of the pipeline
with open("results/cds.txt", "w") as f:
    f.write(f"The HCMV genome (GCF_000845245.1) has {total_cds} CDS.\n")
