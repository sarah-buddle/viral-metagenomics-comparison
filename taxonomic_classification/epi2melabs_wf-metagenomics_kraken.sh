## Run ONT's EPI2ME-labes wf-metagenomics pipeline ##

for sample in ${samples[@]}; do

    # Run analysis
    nextflow run epi2me-labs/wf-metagenomics -profile singularity -w ${results}/work \
    --fastq ${data}/${run}/${sample}_${run}.fq.gz \
    --classifier kraken2 \
    --database $kraken_db \
    --taxonomy $taxdump \
    --out_dir ${results}/epi2me_kraken/${run}/${sample}

    # Combine reports
    samplerun=$(echo -e $sample"\t"$run"\t"$db_name)

    report=${results}/epi2me_kraken/${run}/${sample}/kraken/${sample}_${run}.kreport.txt

    awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report |
    tail -n +2 \
    >> ${results}/epi2me_kraken/all_epi2me_kraken.txt

    report=${results}/epi2me_kraken/${run}/${sample}/bracken/${sample}_${run}.kreport_bracken_species.txt

    awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report |
    tail -n +2 \
    >> ${results}/epi2me_kraken/all_epi2me_bracken.txt

done