## Randomly subsample 5 Gb per sample for Illumina data ##

for sample in ${samples[@]}; do

    $seqtk sample -s100 ${data}/${run}/${sample}_${run}_1.fq.gz 16666667 \
    > ${data}/${run}/${sample}_sub_${run}_1.fq

     $seqtk sample -s100 ${data}/${run}/${sample}_${run}_2.fq.gz 16666667 \
    > ${data}/${run}/${sample}_sub_${run}_2.fq

    gzip ${data}/${run}/${sample}_sub_${run}_1.fq
    gzip ${data}/${run}/${sample}_sub_${run}_2.fq

done





