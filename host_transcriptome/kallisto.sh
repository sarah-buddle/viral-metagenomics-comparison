## Kalisto for human transcript quantification ##

mamba activate kallisto

# Build index
kallisto index -i $human_transcripts ${human_transcripts}.fa

# Quantify transcripts

for sample in ${samples[@]}; do

    mkdir -p ${results}/kallisto/${run}/${sample}

    # Illumina
    kallisto quant -i $human_transcripts \
    -o ${results}/kallisto/${run}/${sample} \
    ${data}/${run}/${sample}_${run}_1.fq.gz \
    ${data}/${run}/${sample}_${run}_2.fq.gz

    # Nanopore
    mean=$(awk '{sum+=$1} END {print sum/NR}' ${results}/lengths/${run}/${sample}_${run}_lengths.txt)
    sd=$(awk '{sum+=$1; sumsq+=$1*$1} END {print sqrt(sumsq/NR - (sum/NR)^2)}' ${results}/lengths/${run}/${sample}_${run}_lengths.txt)

    kallisto quant -i $human_transcripts \
    -o ${results}/kallisto/${run}/${sample} \
    --single ${data}/${run}/${sample}_${run}.fq.gz \
    --fragment-length=$mean \
    --sd=$sd

    awk -F '\t' '{split($1, arr, "|"); printf "%s", arr[1]; for (i = 2; i <= NF; i++) printf "\t%s", $i; printf "\n"}' \
    ${results}/kallisto/${run}/${sample}/abundance.tsv \
    > ${results}/kallisto/${run}/${sample}/abundance_newname.tsv

done
