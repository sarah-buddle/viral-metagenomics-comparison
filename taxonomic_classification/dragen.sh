## Upload data to Basespace and run Illumina's Dragen metagenomics workflow ##

bs auth

# Upload database to databases project
db_project_id=$(bs get project -n Databases --terse)
tar -czvf ${kraken_db}.tar.gz ${kraken_db}
bs upload dataset -p ${db_project_id} --type common.files --name ${kraken_db_name} ${kraken_db}.tar.gz

# Create project
bs create project --name="${run}"

# Upload samples and start analysis
for sample in ${samples[@]}; do

    mkdir -p ${links}/${run}/${sample}

    # Create biosample
    bs create biosample --name="${sample}_${run}" --project="${run}"

    # Find project/sample IDs from names
    project_id=$(bs get project -n $run --terse)
    sample_id=$(bs get biosample -n ${sample}_${run} --terse)
    db_project_id=$(bs get project -n Databases --terse)
    db_file_id=$(bs contents project --id ${db_project_id} | grep "${kraken_db_name}" | awk '{print $2}')

    # Create symlinks to conform to naming conventions
    ln -s ${data}/${run}/${sample}_${run}_1.fq.gz  ${links}/${run}/${sample}/sample_S1_L001_R1_001.fastq.gz
    ln -s ${data}/${run}/${sample}_${run}_2.fq.gz  ${links}/${run}/${sample}/sample_S1_L001_R2_001.fastq.gz

    # Uplaod datasets to basespace
    bs upload dataset --biosample-name="${sample}_${run}" --project=${project_id} \
    ${links}/${run}/${sample}/sample_S1_L001_R1_001.fastq.gz \
    ${links}/${run}/${sample}/sample_S1_L001_R2_001.fastq.gz

    # Launch metagenomics pipeline with custom database
    bs launch application \
    -n "DRAGEN Metagenomics Pipeline" \
    --app-version 3.5.12 \
    -o project-id:${project_id} \
    -o sample-id:${sample_id} \
    -o db:custom \
    -o custom-db:${db_file_id} \
    -o detection-report:none \
    -o basespace-labs-disclaimer:Accepted

done 

# Download results once analysis complete

mkdir ${results}/dragen/${run}/

for sample in ${samples[@]}; do

    # Download
    bs download biosample --name="${sample}_${run}" --extension="tsv" --no-metadata \
    --output ${results}/dragen/${run}/

    mv ${results}/dragen/${run}/${sample}_${run}*/${sample}_${run}.microbe-classification-report.tsv \
    ${results}/dragen/${run}/

    # Combine reports
    samplerun=$(echo -e $sample"\t"$run"\t"$kraken_db_name)

    report=${results}/dragen/${run}/${sample}_${run}.microbe-classification-report.tsv

    awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report |
    tail -n +2 \
    >> ${results}/dragen/all_dragen.txt

done

 