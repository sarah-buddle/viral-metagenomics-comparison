## Run MEGAN-LR through pb-metagenomics pipeline ##

## Please see: https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/master/docs/Tutorial-Taxonomic-Profiling-Minimap-Megan.md

#### Setup ####

# Remove previous sample names from config files
grep -v '^  -' ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan/configs/Sample-Config.yaml \
> ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan/configs/Sample-Config-temp.yaml

mv ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan/configs/Sample-Config-temp.yaml \
${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan/configs/Sample-Config.yaml

grep -v '^  -' ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan/configs/Sample-Config.yaml \
> ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan/configs/Sample-Config-temp.yaml

mv ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan/configs/Sample-Config-temp.yaml \
${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan/configs/Sample-Config.yaml

for sample in ${samples[@]}; do

    # Convert fastq (output of taxprofiler preprocessing and host removal pipeline) to fasta for input
    zcat ${results}/samtools/fastq/${sample}_${run}.unmapped_other.fastq.gz |
    paste - - - - | cut -f 1,2 | sed 's/^@/>/'  | tr "\t" "\n" \
    > ${results}/samtools/fastq/${sample}_${run}.fasta

    # Create a symlink to the input file in the inputs directory
    ln -s ${results}/samtools/fastq/${sample}_${run}.fasta \
    ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan/inputs/${sample}_${run}.fasta

    ln -s ${results}/samtools/fastq/${sample}_${run}.fasta \
    ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan/inputs/${sample}_${run}.fasta

    # Add sample names to config files
    echo -e "  - \""${sample}_${run}"\"" \
    >> ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan/configs/Sample-Config.yaml

    echo -e "  - \""${sample}_${run}"\"" \
    >> ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan/configs/Sample-Config.yaml

done

### Nucleotide mode ####

# Run analysis with snakemake 
conda activate snakemake
cd ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan

snakemake --snakefile ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan/Snakefile-minimap-megan-v2 \
--configfile ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Minimap-Megan/configs/Sample-Config.yaml \
--use-conda --cores 8

date

#### Protein mode ####

# Make diamond database for protein mode
cd /SAN/breuerlab/BTRU-scratch2/sarah/software/diamond_refseq_db
mamba activate diamond

conda activate snakemake
cd ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan

snakemake --snakefile ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan/Snakefile-diamond-megan \
--configfile ${software}/pb-metagenomics-tools/Taxonomic-Profiling-Diamond-Megan/configs/Sample-Config.yaml \
--use-conda --cores 8

# Combine reports once analysis complete
for sample in ${samples[@]}; do

    db_name=refseq-2023-06-08-nucleotide-v2

    report=${megan_dir_nt}/${sample}_${run}.diamond_megan.kreport.unfiltered.txt

    samplerun=$(echo -e $sample "\t" $run "\t" $db_name)

    awk -v n="$samplerun"  'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report \
    >> ${results}/megan/all_megan.txt

    db_name=refseq-2023-06-08-protein-v2

    report=${megan_dir_pr}/${sample}_${run}.diamond_megan.kreport.unfiltered.txt

    samplerun=$(echo -e $sample "\t" $run "\t" $db_name)

    awk -v n="$samplerun"  'BEGIN{FS=OFS="\t"}{print n OFS $0}' $report \
    >> ${results}/megan/all_megan.txt

done

