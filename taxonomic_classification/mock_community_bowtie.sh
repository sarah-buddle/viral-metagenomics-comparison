## Align to reference genomes of viruses in mock community ##

declare -A taxons=( ["Human mastadenovirus F"]=130309 ["Human herpesvirus 5"]=10359 ["Human respiratory syncytial virus"]=11250 \
 ["Influenza B virus"]=11520 ["Reovirus 3"]=351073 \
 ["Zika virus"]=64320 ["Lambda phage"]=2681611 ["MS2 phage"]=12022 )

# Build bowtie indexes
for i in "${!taxons[@]}"; do

    underscore_organism=$(echo "$organism" | tr ' ' '_')

    ${bowtie}-build ${genome_dir}/${community}/Genomes/${underscore_organism}_genome.fasta \
    ${genome_dir}/${community}/Genomes/${underscore_organism}_genome

done

community=msa_2008

mkdir ${results}/bowtie/${run}

for sample in ${samples[@]}; do

    for i in "${!taxons[@]}"; do

        organism=$i
        taxon_id=${taxons[$i]}
        underscore_organism=$(echo "$organism" | tr ' ' '_')

        # Align with bowtie
        (${bowtie} -x ${genome_dir}/${community}/Genomes/${underscore_organism}_genome \
        --very-sensitive \
        -1 ${results}/bowtie2/align/${sample}_${run}.unmapped_1.fastq.gz \
        -2 ${results}/bowtie2/align/${sample}_${run}.unmapped_2.fastq.gz \
        -S ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.sam) \
        2> ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_log.txt

        # Extract mapped reads
        ${samtools} view -bF 4 -h ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.sam |
        ${samtools} sort - \
        > ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.bam

        rm ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.sam

        # Remove duplicates
        $samtools collate -@ 4 -O -u ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.bam |
        $samtools fixmate -@ 4 -m -u - - |
        $samtools sort -@ 4 -u - |
        $samtools markdup -@ 4 -r -d 2500 - ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_dedup.bam \
        -f ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_dedup_stats.txt 

        # Index
        ${samtools} index ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_dedup.bam

        # Depth (deduplicated data)
        ${samtools} depth ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_dedup.bam \
        > ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_depth.txt

        # Depth (non-deduplicated data)
        ${samtools} depth ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.bam \
        > ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_nodedup_depth.txt

        # Extract fastq (deduplicated data)
        ${samtools} fastq ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_dedup.bam \
        > ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.fastq

        # Extract fastq (deduplicated data)
        ${samtools} fastq ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.bam \
        > ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_nodedup.fastq

        # Make report - read counts
        bowtie_reads=$(awk 'NR % 4 == 1' ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.fastq |
        wc -l | awk '{print $1}')

        bowtie_bp=$(awk 'NR % 4 == 2' ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.fastq | 
        tr -d '\n' | tr -d '[:space:]' | wc -m)

        echo -e "bowtie\t"$run"\t"$sample"\t"$organism"\t"$bowtie_reads"\t"$bowtie_bp \
        >> ${results}/bowtie/all_bowtie.txt

        bowtie_reads_nodedup=$(awk 'NR % 4 == 1' ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_nodedup.fastq |
        wc -l | awk '{print $1}')

        bowtie_bp_nodedup=$(awk 'NR % 4 == 2' ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_nodedup.fastq | 
        tr -d '\n' | tr -d '[:space:]' | wc -m)

        echo -e "bowtie\t"$run"\t"$sample"\t"$organism"\t"$bowtie_reads_nodedup"\t"$bowtie_bp_nodedup \
        >> ${results}/bowtie/all_bowtie_nodedup.txt

        # Coverage (deduplicated data)
        ${samtools} coverage ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_dedup.bam \
        > ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_coverage.txt

        # Make coverage report
        samplerun=$(echo -e $sample "\t" $run "\t" $organism)

        awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' \
        ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_coverage.txt |
        tail -n +2 \
        >> ${results}/bowtie/coverage.txt

        ${samtools} coverage ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}.bam \
        > ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_nodedup_coverage.txt

        samplerun=$(echo -e $sample "\t" $run "\t" $organism)

        awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' \
        ${results}/bowtie/${run}/${sample}_${run}_${underscore_organism}_nodedup_coverage.txt |
        tail -n +2 \
        >> ${results}/bowtie/coverage_nodedup.txt

        done

    done

done