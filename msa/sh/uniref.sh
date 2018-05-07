mkdir tmp
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz
gunzip uniref50.fasta.gz
mmseqs createdb uniref50.fasta uniref100.db
mmseqs createindex uniref50.db tmp
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz
mmseqs createdb uniref90.fasta uniref100.db
mmseqs createindex uniref90.db tmp
curl -O ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
gunzip uniref100.fasta.gz
mmseqs createdb uniref100.fasta uniref100.db
mmseqs createindex uniref100.db tmp
rm -rf tmp
