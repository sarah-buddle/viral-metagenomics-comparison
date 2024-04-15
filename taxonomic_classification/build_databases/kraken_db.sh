## Build Kraken2 database ##

mamba activate kraken

# Add taxdump and accession2taxid files to taxonomy dir first

kraken2-build --add-to-library $fasta --db refseq-2023-06-08-nucleotide-v2

kraken2-build --build --db refseq-2023-06-08-nucleotide-v2 --threads 4

bracken-build -d refseq-2023-06-08-nucleotide-v2 -t 4

