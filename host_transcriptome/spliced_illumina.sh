## Extract reads mapping across splice junctions and quantify transcirpts with kallisto ##

mamba activate kallisto

$star runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $index \
--genomeFastaFiles /SAN/breuerlab/BTRU-scratch/sarah/data/human_genome/ensembl_GRCh38_107/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /SAN/breuerlab/BTRU-scratch/sarah/data/human_genome/ensembl_GRCh38_107/gtf/Homo_sapiens.GRCh38.107.gtf \
--sjdbOverhang 149

for sample in ${samples[@]}; do

    # QC
    mkdir ${results}/fastp_manual/${run}

    $fastp -i ${data}/${run}/${sample}_${run}_1.fq.gz -I ${data}/${run}/${sample}_${run}_2.fq.gz \
    -o ${results}/fastp_manual/${run}/${sample}_${run}_1.fq.gz -O ${results}/fastp_manual/${run}/${sample}_${run}_2.fq.gz \
    -h ${results}/fastp_manual/${run}/${sample}_${run}_fastp.html -j ${results}/fastp_manual/${run}/${sample}_${run}_fastp.json

    # Align with STAR

    mkdir -p ${results}/star/${run}/${sample}

    $star runThreadN 8 \
    --genomeDir $index \
    --readFilesIn ${results}/fastp_manual/${run}/${sample}_${run}_1.fq.gz ${results}/fastp_manual/${run}/${sample}_${run}_2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix ${results}/star/${run}/${sample}/${sample}_${run}_

    $samtools view -bh ${results}/star/${run}/${sample}/${sample}_${run}_Aligned.out.sam |
    $samtools sort - \
    > ${results}/star/${run}/${sample}/${sample}_${run}_Aligned.out.bam

    $samtools index ${results}/star/${run}/${sample}/${sample}_${run}_Aligned.out.bam

    # This looks for the cigar string (N - how sam/bam represents gapped alignments) in the alignment
    awk '$0 ~ /^@/ || $6 ~ /N/' ${results}/star/${run}/${sample}/${sample}_${run}_Aligned.out.sam |
    $samtools view -bh > ${results}/star/${run}/${sample}/${sample}_${run}_spliced.bam

    # Convert to fastq and extract the pairs of singleton reads
    $samtools fastq ${results}/star/${run}/${sample}/${sample}_${run}_spliced.bam \
    -1 ${results}/star/${run}/${sample}/${sample}_${run}_spliced_p1.fastq \
    -2 ${results}/star/${run}/${sample}/${sample}_${run}_spliced_p2.fastq \
    -s ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s.fastq

    grep '@' ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s.fastq |
    sed 's/@//' \
    > ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s_ids.txt

    $seqtk subseq ${results}/fastp_manual/${run}/${sample}_${run}_1.fq.gz \
    ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s_ids.txt \
    > ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s1.fastq

    $seqtk subseq ${results}/fastp_manual/${run}/${sample}_${run}_2.fq.gz \
    ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s_ids.txt \
    > ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s2.fastq

    cat ${results}/star/${run}/${sample}/${sample}_${run}_spliced_p1.fastq \
    ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s1.fastq \
    > ${results}/star/${run}/${sample}/${sample}_${run}_spliced_1.fastq

    cat ${results}/star/${run}/${sample}/${sample}_${run}_spliced_p2.fastq \
    ${results}/star/${run}/${sample}/${sample}_${run}_spliced_s2.fastq \
    > ${results}/star/${run}/${sample}/${sample}_${run}_spliced_2.fastq

    gzip -f ${results}/star/${run}/${sample}/*.fastq
    rm ${results}/star/${run}/${sample}/${sample}_${run}_Aligned.out.sam

    # Quantify transcripts with kallisto
    mkdir -p ${results}/kallisto_spliced/${run}/${sample}

    kallisto quant -i $human_transcripts \
    -o ${results}/kallisto_spliced/${run}/${sample} \
    ${results}/star/${run}/${sample}/${sample}_${run}_spliced_1.fastq.gz \
    ${results}/star/${run}/${sample}/${sample}_${run}_spliced_2.fastq.gz

    awk -F '\t' '{split($1, arr, "|"); printf "%s", arr[1]; for (i = 2; i <= NF; i++) printf "\t%s", $i; printf "\n"}' \
    ${results}/kallisto_spliced/${run}/${sample}/abundance.tsv \
    > ${results}/kallisto_spliced/${run}/${sample}/abundance_newname.tsv

done

