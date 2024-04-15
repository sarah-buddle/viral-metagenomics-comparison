## Randomly subsample 5 Gb per sample for Nanopore data ##

# Nanopore
for sample in ${samples[@]}; do

    zcat ${data}/${run}/${sample}_${run}.fq.gz |s
    paste - - - - | shuf | paste -d'\t' - - - - | tr '\t' '\n' |
    head -n -2 \
    > ${data}/${run}/${sample}_${run}_shuffled.fq

    python3 ${code}/subsample_fastq.py \
    ${data}/${run}/${sample}_${run}_shuffled.fq \
    ${data}/${run}/${sample}_sub_${run}.fq \
    5000000000

    cat ${data}/${run}/${sample}_sub_${run}.fq |
    awk  'NR%4==2 {print length}' | 
    sort -nr \
    > ${results}/lengths/${run}/${sample}_sub_${run}_lengths.txt

    N_READS=$(wc -l ${results}/lengths/${run}/${sample}_sub_${run}_lengths.txt | awk '{print $1}')

    TOTAL_LENGTH=$(awk '{SUM+=$1}END{print SUM}' ${results}/lengths/${run}/${sample}_sub_${run}_lengths.txt)

    rm ${data}/${run}/${sample}_${run}_shuffled.fq
    gzip ${data}/${run}/${sample}_sub_${run}.fq

done

