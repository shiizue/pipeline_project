from Bio import Entrez, SeqIO

# CHANGE THIS WITH CONFIG FILE LATER*
Entrez.email = "skyel1005@gmail.com"
refseq = "NC_006273.2"

handle = Entrez.efetch(db="nucleotide", id=refseq, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")

cds_list = []
total_cds = 0


for feature in record.features:
    if feature.type == "CDS":
        total_cds+=1
        protein_id=feature.qualifiers["protein_id"][0]
        cds_seq=feature.extract(record.seq)
        cds_list.append([protein_id,cds_seq])

with open("hcmv_cds.fasta", "w") as f:
    for protein_id, cds_seq in cds_list:
        f.write(f">{protein_id}\n{cds_seq}\n")

with open("PipelineReport.txt", "w") as f:
    f.write(f"The HCMV genome (GCF_000845245.1) has {total_cds} CDS.\n")