### Preprocess Nanopore data ###

for sample in ${samples[@]}; do

    zcat ${data}/${run}/fastq_pass/${sample}/*.fastq.gz > ${data}/${run}/${sample}_${run}.fq

    gzip ${data}/${run}/${sample}_${run}.fq

done

