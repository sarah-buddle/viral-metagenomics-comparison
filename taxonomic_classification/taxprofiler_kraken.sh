## Run Kraken2 + Bracken through Taxprofiler ##

# Run taxprofiler - use preprocessed reads as input
nextflow run nf-core/taxprofiler -r dev -profile singularity -with-tower -w ${results}/work --max_cpus 1 \
--input ${results}/samplesheets/samplesheet_filtered_${run}.csv \
--databases ${results}/databases/database_kraken.csv \
--outdir ${results} \
--run_kraken2 --kraken2_save_readclassifications \
--run_bracken

# For Nanopore - run bracken separately since not done through taxprofiler
for sample in ${samples[@]}; do

    bracken -d $db \
    -i ${results}/kraken2/${db_name}/${sample}_se_${run}_${db_name}.kraken2.kraken2.report.txt \
    -o ${results}/bracken/${db_name}/${sample}_se_${run}_${db_name}.bracken.tsv
    
    # Create report
    report=${results}/kraken2/${db_name}/${sample}_se_${run}_${db_name}.kraken2.kraken2.report.txt

    samplerun=$(echo -e $sample "\t" $run "\t" $db_name)

    awk -v n="$samplerun"  'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report \
    >> ${results}/kraken2/all_kraken.txt

    report=${results}/bracken/${db_name}/${sample}_se_${run}_${db_name}.bracken.tsv

    awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report |
    tail -n +2 \
    >> ${results}/bracken/all_bracken.txt

done

# Create report - Illumina
for sample in ${samples[@]}; do
 
    report=${results}/kraken2/${db_name}/${sample}_pe_${run}_${db_name}.kraken2.kraken2.report.txt

    samplerun=$(echo -e $sample "\t" $run "\t" $db_name)

    awk -v n="$samplerun"  'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report \
    >> ${results}/kraken2/all_kraken.txt

    report=${results}/bracken/${db_name}/${sample}_pe_${run}_${db_name}.bracken.tsv

    awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report |
    tail -n +2 \
    >> ${results}/bracken/all_bracken.txt

done
