## Timings plot ##

library(tidyverse)

# Import start times and read lengths from sequencing summary file
timings_raw_dna <- read.delim(input_dna, sep = "\t", header = TRUE) %>% 
  mutate(platform = "nanopore_dna", output = NA)

timings_raw_rna1 <- read.delim(input_rna1, sep = "\t", header = TRUE) %>% 
  mutate(platform = "nanopore_rna", output = NA)

# RNA run was interrupted so had two separate seq sum files
max_rna1 <- max(timings_raw_rna1$start_time)

timings_raw_rna2 <- read.delim(input_rna2, sep = "\t", header = TRUE) %>% 
  mutate(platform = "nanopore_rna", output = NA)
  
timings_rna2 <- timings_raw_rna2 %>% 
  mutate(start_time = start_time + max_rna1)

# Illumina kits data volumes and timins
illumina <- read.csv(illumina_timings) %>%
  select(-c(cycles, reads_passing_filter)) %>% 
  mutate(total_bases_passfail = NA)

# Calculate cumulative data volumes from Nanopore
timings <- rbind(timings_raw_dna, timings_raw_rna1, timings_rna2) %>% 
  dplyr::filter(passes_filtering == TRUE) %>% 
  dplyr::rename(read_length = sequence_length_template) %>% 
  dplyr::group_by(platform) %>% 
  dplyr::mutate(total_bases = cumsum(as.numeric(read_length)),
        time_hours = start_time / 3600) %>% 
  dplyr::group_by(platform, passes_filtering) %>% 
  dplyr::mutate(total_bases_passfail = cumsum(as.numeric(read_length))) %>% 
  dplyr::filter(row_number() %% 1000 == 0) %>% 
  dplyr::select(-c(start_time, read_length)) %>% 
  dplyr::mutate(total_gb = NA) %>% 
  rbind(illumina)

# Plot
timings %>% 
  dplyr::filter(time_hours <= 73) %>% 
  ggplot2::ggplot() +
  geom_smooth(aes(x = time_hours, y = total_bases_passfail/10^9), color = "blue") +
  geom_point(aes(x = time_hours, y = total_gb), size = 2.5, colour = "darkturquoise") +
  theme(panel.background = element_blank(),
       panel.border = element_rect(colour = "black", fill = NA),
       legend.position = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = seq(0, 72, 12)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  ylab("Total Gb") +
  xlab("Time (hours)")

ggsave(output, height = 3, width = 4, dpi = 750)
