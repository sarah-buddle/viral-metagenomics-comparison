## Run BLAST and DIAMOND on taxprofiler preprocessing/host removal output before input to metaMix ##

# ONT method

for sample in ${samples[@]}; do

    echo $sample
    date

    # Convert to fasta
    zcat ${results}/samtools/fastq/${sample}_${run}.unmapped_other.fastq.gz |
    $seqtk seq -a \
    > ${results}/samtools/fastq/${sample}_${run}.fasta

    mkdir -p ${results}/metamix/${run}/${sample}_${run}/nucleotide
    mkdir ${results}/metamix/${run}/${sample}_${run}/metamix_nucleotide_norRNA

    # Run blastn
    $blastn -db $dbnucl -query ${results}/samtools/fastq/${sample}_${run}.fasta \
    -outfmt 6 -max_target_seqs 10 -max_hsps 1 -num_threads 6 \
    > ${results}/metamix/${run}/${sample}_${run}/nucleotide/${sample}_${run}_megaBLAST_norRNA.tab

    # Create read lengths file
    awk '/^>/ {if (NR>1) print ""; printf "%s\t", substr($0,2); next} {printf "%s", $0} END {print ""}' \
    ${results}/samtools/fastq/${sample}_${run}.fasta | 
    awk -F '\t' '{print $1 "\t" length($2)}' \
    > ${results}/metamix/${run}/${sample}_${run}/metamix_nucleotide_norRNA/read_lengths.tab

    mkdir ${results}/metamix/${run}/${sample}_${run}/protein
    mkdir ${results}/metamix/${run}/${sample}_${run}/metamix_protein

    cp ${results}/metamix/${run}/${sample}_${run}/metamix_nucleotide_norRNA/read_lengths.tab \
    ${results}/metamix/${run}/${sample}_${run}/metamix_protein/read_lengths.tab

    # Run diamond
    $diamond blastx -d $dbprot \
    -q ${results}/samtools/fastq/${sample}_${run}.fasta \
    -f 6  -k 10 -p 6 -sensitive -e 1 \
    > ${results}/metamix/${run}/${sample}_${run}/protein/${sample}_${run}_diamond.tab

done
