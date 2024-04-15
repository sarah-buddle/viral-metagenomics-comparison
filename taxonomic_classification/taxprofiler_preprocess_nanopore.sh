## QC and host filtering with Taxprofiler ##

nextflow run nf-core/taxprofiler -r dev -profile conda -with-tower -w ${results}/work --max_cpus 4 \
-c ${software}/nf-core-configs/custom_resources.conf \
--input ${results}/samplesheets/samplesheet.csv \
--databases ${results}/databases/database_kraken_refseq_nucleotide_v2.csv \
--outdir ${results} \
--perform_longread_qc \
--longread_qc_skipqualityfilter \
--perform_longread_hostremoval \
--hostremoval_reference ${human_genome} \
--save_hostremoval_index \
--save_hostremoval_unmapped