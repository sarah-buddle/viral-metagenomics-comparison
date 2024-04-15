# Make diamond database for MEGAN-LR and metaMix protein mode
mamba activate diamond

diamond makedb --in ${dbcustomprot} --db refseq-2023-06-08-protein --threads 8