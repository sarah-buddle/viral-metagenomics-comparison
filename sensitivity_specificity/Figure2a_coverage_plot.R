### Coverage ####

library(tidyverse)

#### Levels and labels ####

dilutions <- c(10, 100, 1000, 10000)

run_levels <- c("twist_101123", "illumina_250923", "nanopore_270923")

run_labels <- c("Twist Viral Research Panel (VRP)", "Untargeted Illumina", "Untargeted ONT")

run_colours <-setNames(c("#619CFF", "#F8766D", "#00BA38"),
                       c("Twist Viral Research Panel (VRP)", "Untargeted Illumina", "Untargeted ONT"))


format_labels <- function(l) {
  
  label <- 6 * 10^l
  
  parse(text=label)
  
}


msa_2008_levels <- c("Human betaherpesvirus 5", "Human herpesvirus 5", "Cytomegalovirus humanbeta5",
                     "Human mastadenovirus F", 
                     "Human orthopneumovirus", "Orthopneumovirus hominis", "Human respiratory syncytial virus",
                     "Influenza B virus", "Betainfluenzavirus influenzae",
                     "Mammalian orthoreovirus", "Reovirus 3",
                     "Zika virus", "Orthoflavivirus zikaense",
                     "Lambda phage", "Escherichia virus Lambda", "Lambdavirus lambda",
                     "MS2 phage", "Escherichia virus MS2", "Emesvirus zinderi")

msa_2008_labels_oneline <- c("Human betaherpesvirus 5", "Human betaherpesvirus 5", "Human betaherpesvirus 5",
                             "Human mastadenovirus F", 
                             "Human orthopneumovirus", "Human orthopneumovirus", "Human orthopneumovirus",
                             "Influenza B virus", "Influenza B virus",
                             "Mammalian orthoreovirus", "Mammalian orthoreovirus", 
                             "Zika virus", "Zika virus",
                             "Lambda phage", "Lambda phage", "Lambda phage",
                             "MS2 phage", "MS2 phage", "MS2 phage")

msa_2008_labels <- c("Human\nbetaherpesvirus 5", "Human\nbetaherpesvirus 5", "Human\nbetaherpesvirus 5",
                     "Human\nmastadenovirus F", 
                     "Human\northopneumovirus", "Human\northopneumovirus", "Human\northopneumovirus",
                     "Influenza B\nvirus", "Influenza B\nvirus",
                     "Mammalian\northoreovirus", "Mammalian\northoreovirus", 
                     "Zika\nvirus", "Zika\nvirus",
                     "Lambda\nphage", "Lambda\nphage", "Lambda\nphage",
                     "MS2\nphage", "MS2\nphage", "MS2\nphage")

#### File inputs ####
controls <- read.csv(controls)

community_taxid_mapping <- read.csv(community_taxid_mapping)

genome_lengths <- read.csv(genome_lengths) %>% 
  dplyr::mutate(species = factor(species, msa_2008_levels, msa_2008_labels)) %>% 
  dplyr::rename(name = species)

#### Coverage plot ####

# Import data
coverage_minimap <- read.delim(input_coverage_minimap, sep = "\t", header = TRUE, row.names = NULL) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(length = endpos - startpos)

coverage_bowtie <- read.delim(input_coverage_bowtie, sep = "\t", header = TRUE, row.names = NULL) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(length = endpos - startpos)

# Clean
average_coverage <- rbind(coverage_bowtie, coverage_minimap) %>% 
  dplyr::filter(grepl("_sub", sample) & run %in% run_levels) %>% 
  dplyr::select(-c(startpos, endpos, numreads, covbases)) %>% 
  dplyr::group_by(sample, run, name) %>% 
  dplyr::summarise(coverage = weighted.mean(coverage, length), depth = weighted.mean(meandepth, length), 
            baseq = weighted.mean(meanbaseq, length), mapq = weighted.mean(meanmapq, length)) %>% 
  dplyr::arrange(run, sample, name) %>% 
  dplyr::mutate(name = gsub("_", " ", name)) %>% 
  dplyr::left_join(controls, by = c("run", "sample")) %>% 
  dplyr::select(-control) %>% 
  dplyr::mutate(dilution = as.character(dilution)) %>% 
  dplyr::mutate(name = factor(name, msa_2008_levels, msa_2008_labels),
         run = factor(run, run_levels, run_labels)) %>% 
  dplyr::filter((name %in% c("Human\nbetaherpesvirus 5", "Human\nmastadenovirus F") & dnarna == "DNA" & run != "Twist Viral Research Panel (VRP)") |
           (name %in% c("Human\northopneumovirus", "Influenza B\nvirus", "Mammalian\northoreovirus",  "Zika\nvirus") 
            & dnarna == "RNA" & run != "Twist Viral Research Panel (VRP)") |
           (!(name %in% c("Lambda\nphage", "MS2\nphage")) & dnarna == "DNARNA" & run == "Twist Viral Research Panel (VRP)") &
           !is.na(dilution)) %>% 
  dplyr::left_join(genome_lengths, by = "name") %>% 
  dplyr::mutate(depth_norm = depth / genome_length)

# Line plot coverage
average_coverage %>% 
  dplyr::filter(dilution %in% dilutions) %>% 
  dplyr::group_by(run, name, repeat_group, dilution) %>% 
  dplyr::summarise(mean_coverage = mean(coverage), min_coverage = min(coverage), max_coverage = max(coverage)) %>% 
  dplyr::mutate(dilution = as.numeric(dilution)) %>% 
  ggplot2::ggplot(mapping = aes(x = log10(dilution), y = mean_coverage, color = run)) +
  geom_smooth(lwd = 0.8) +
  facet_wrap(vars(name))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0)) +
  scale_x_continuous(breaks = c(1,2,3,4), labels = format_labels) +
  scale_color_manual(values = run_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0, 100)) +
  geom_errorbar(aes(ymin = min_coverage, max = max_coverage), width = 0.1) +
  ylab("Percentage genome coverage") +
  xlab("Genome copies per ml (gc/ml)")

ggsave(paste0(outdir, "/coverage.png"), height = 3.5, width = 6)



