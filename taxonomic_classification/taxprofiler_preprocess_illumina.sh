## QC and host filtering with Taxprofiler ##

nextflow run nf-core/taxprofiler -r dev -profile singularity -with-tower -w ${results}/work --max_cpus 4 \
--input ${results}/samplesheets/samplesheet_${run}.csv --outdir ${results} \
--databases ${results}/databases/database_kraken.csv \
--perform_shortread_qc \
--perform_shortread_hostremoval \
--hostremoval_reference ${human_genome} \
--shortread_hostremoval_index ${human_index_dir_bowtie} \
--save_hostremoval_unmapped

