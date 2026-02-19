import_os

download_data = (
    "datasets download virus genome taxon Betaherpesvirinae --refseq --include genome"
)
unzip_data = "unzip ncbi_dataset.zip"
make_blast_db = "makeblastdb -in betaherpesvirinae.fasta -out blast/betaherpesvirinae -title betaherpesvirina -dbtype nucl"

os.system(download_data)
os.system(unzip_data)
os.system(make_blast_db)
