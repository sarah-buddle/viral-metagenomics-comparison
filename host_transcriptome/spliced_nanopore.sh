## Extract reads mapping across splice junctions and quantify transcirpts with kallisto ##

mamba activate kallisto

$minimap -d ${human_genome}.mmi ${human_genome}.fa

for sample in ${samples[@]}; do

    # Align with minimap
    mkdir -p ${results}/minimap_human/${run}/${sample}

    $minimap -ax map-ont --splice ${human_genome}.fa \
    ${data}/${run}/${sample}_${run}.fq.gz \
    > ${results}/minimap_human/${run}/${sample}/${sample}_${run}_Aligned.out.sam

    $samtools view -bh ${results}/minimap_human/${run}/${sample}/${sample}_${run}_Aligned.out.sam |
    $samtools sort - \
    > ${results}/minimap_human/${run}/${sample}/${sample}_${run}_Aligned.out.bam

    $samtools index ${results}/minimap_human/${run}/${sample}/${sample}_${run}_Aligned.out.bam

    # This looks for the cigar string (N - how sam/bam represents gapped alignments) in the alignment
    awk '$0 ~ /^@/ || $6 ~ /N/' ${results}/minimap_human/${run}/${sample}/${sample}_${run}_Aligned.out.sam |
    $samtools view -bh > ${results}/minimap_human/${run}/${sample}/${sample}_${run}_spliced.bam

    # Convert to fastq and extract the pairs of singleton reads
    $samtools fastq ${results}/minimap_human/${run}/${sample}/${sample}_${run}_spliced.bam \
    > ${results}/minimap_human/${run}/${sample}/${sample}_${run}_spliced.fastq

    gzip -f ${results}/minimap_human/${run}/${sample}/${sample}_${run}_spliced.fastq
    rm ${results}/minimap_human/${run}/${sample}/${sample}_${run}_Aligned.out.sam

    # Quantify transcripts with kallisto
    mkdir -p ${results}/kallisto_spliced/${run}/${sample}

    mean=$(awk '{sum+=$1} END {print sum/NR}' ${results}/lengths/${run}/${sample}_${run}_lengths.txt)
    sd=$(awk '{sum+=$1; sumsq+=$1*$1} END {print sqrt(sumsq/NR - (sum/NR)^2)}' ${results}/lengths/${run}/${sample}_${run}_lengths.txt)

    kallisto quant -i $human_transcripts \
    -o ${results}/kallisto_spliced/${run}/${sample} \
    --single ${results}/minimap_human/${run}/${sample}/${sample}_${run}_spliced.fastq.gz \
    --fragment-length=$mean \
    --sd=$sd

    awk -F '\t' '{split($1, arr, "|"); printf "%s", arr[1]; for (i = 2; i <= NF; i++) printf "\t%s", $i; printf "\n"}' \
    ${results}/kallisto_spliced/${run}/${sample}/abundance.tsv \
    > ${results}/kallisto_spliced/${run}/${sample}/abundance_newname.tsv

done

