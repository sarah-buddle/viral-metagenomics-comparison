## Align to reference genomes of viruses in mock community ##

declare -A taxons=( ["Human mastadenovirus F"]=130309 ["Human herpesvirus 5"]=10359 ["Human respiratory syncytial virus"]=11250 \
 ["Influenza B virus"]=11520 ["Reovirus 3"]=351073 \
 ["Zika virus"]=64320 ["Lambda phage"]=2681611 ["MS2 phage"]=12022)

# Make minimap2 indexes
for i in "${!taxons[@]}"; do

    underscore_organism=$(echo "$organism" | tr ' ' '_')

    $minimap -d ${genome_dir}/${community}/Genomes/${underscore_organism}_genome.mmi \
    ${genome_dir}/${community}/Genomes/${underscore_organism}_genome.fasta

done

for sample in  ${samples[@]}; do

    for i in "${!taxons[@]}"; do 
        organism=$i
        taxon_id=${taxons[$i]}
        underscore_organism=$(echo "$organism" | tr ' ' '_')

        mkdir -p ${results}/minimap/${run}/${sample}/${underscore_organism}

        # Align to reference genome
        ${minimap} -ax map-ont ${genome_dir}/${community}/Genomes/${underscore_organism}_genome.fasta  \
        ${results}/samtools/fastq/${sample}_${run}.unmapped_other.fastq.gz \
        > ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}.sam

        # Exract aligned reads
        ${samtools} view -bF 4 -h ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}.sam |
        ${samtools} sort \
        > ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.bam

        if test -s ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.bam; then
            rm ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}.sam
        fi

        # Remove duplicates
        $samtools markdup -@ 4 -r -l 500000 ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.bam \
         ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_dedup.bam \
        -f ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_dedup_stats.txt

        # Index bam file
        ${samtools} index ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_dedup.bam 

        # Calculate depth and coverage
        ${samtools} depth ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_dedup.bam \
        > ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_depth.txt

        ${samtools} coverage ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_dedup.bam \
        > ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_dedup_coverage.txt

         ${samtools} depth ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.bam \
        > ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_nodedup_depth.txt

        ${samtools} coverage ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.bam \
        > ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_nodedup_coverage.txt

        # Make fastqs
        ${samtools} fastq ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_dedup.bam  \
        > ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.fastq

        ${samtools} fastq ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.bam  \
        > ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_nodedup.fastq

        samplerun=$(echo -e $sample "\t" $run "\t" $underscore_organism)

        # Make coverage report
        awk -v n="$samplerun" 'BEGIN{FS=OFS="\t"}{print n OFS $0}' \
        ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped_dedup_coverage.txt |
        tail -n +2 \
        >> ${results}/minimap/coverage.txt

        # Make alignment report
        minimap_reads=$(awk 'NR % 4 == 1' ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.fastq |
        wc -l | awk '{print $1}')

        minimap_bp=$(awk 'NR % 4 == 2' ${results}/minimap/${run}/${sample}/${underscore_organism}/${run}_${sample}_${underscore_organism}_mapped.fastq | 
        tr -d '\n' | tr -d '[:space:]' | wc -m)

        echo -e "minimap\t"$run"\t"$sample"\t"$organism"\t"$minimap_reads"\t"$minimap_bp \
        >> ${results}/minimap/all_minimap.txt

        done

done

