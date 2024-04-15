# Extract length and start times from sequencing summary
cut -f 10,16 ${data}/${run}/seq_sum/${id}.txt \
> ${results}/nanopore_timings/${run}/${id}_times.txt