## Run kaiju using taxprofiler ##

# Use preprocessed reads as input
nextflow run nf-core/taxprofiler -r dev -profile singularity -with-tower -w ${results}/work --max_cpus 1 \
--input ${results}/samplesheets/samplesheet_${run}.csv \
--databases ${results}/databases/database_kaiju.csv \
--outdir ${results} \
--run_kaiju \
--kaiju_expand_viruses

# Make report
for sample in ${samples[@]}; do
 
    # Illumina
    report=${results}/kaiju/${db_name}/${sample}_pe_${run}_${db_name}.kaijutable.txt

    # Nanopore
    report=${results}/kaiju/${db_name}/${sample}_se_${run}_${db_name}.kaijutable.txt

    samplerun=$(echo -e $sample "\t" $run "\t" $db_name)c

    awk -v n="$samplerun"  'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report |
    tail -n +2 \
    >> ${results}/kaiju/all_kaiju.txt

done