library(tidyverse)

##### Levels and labels ####

dilutions <- c(10, 100, 1000, 10000)

run_levels <- c("twist_101123", "illumina_250923", "nanopore_270923")

run_labels <- c("Twist Viral Research Panel (VRP)", "Untargeted Illumina", "Untargeted ONT")

run_colours <-setNames(c("#619CFF", "#F8766D", "#00BA38"),
                       c("Twist Viral Research Panel (VRP)", "Untargeted Illumina", "Untargeted ONT"))
                                        


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


tool_levels_full <- c("bowtie_nodedup", "minimap_nodedup", "kraken", "bracken", "dragen", "epi2me_bracken", "megan_nucleotide", "megan_protein", "kaiju_refseq-2023-06-08-nucleotide-v2",
                      "metamix_nucleotide", "metamix_step2_nucleotide", "metamix_protein", "metamix_step2_protein", "czid_nt_count", "czid_nr_count", "onecodex")

tool_labels_full <- c( "Bowtie2/\nMinimap2", "Bowtie2/\nMinimap2",
                      "Kraken2", "Bracken", "Dragen/\nEPI2ME", "Dragen/\nEPI2ME", "MEGAN-LR\n(nucleotide)", "MEGAN-LR\n(protein)", "Kaiju",
                      "metaMix\n(nucleotide)", "metaMix Fast\n(nucleotide)", "metaMix \n(protein)",  "metaMix Fast\n(protein)", 
                      "CZ ID*\n(nucleotide)", "CZ ID*\n(protein)", "One Codex*")

tool_labels_full_oneline <- c( "Bowtie2", "Minimap2",
                       "Kraken2", "Bracken", "Dragen/ EPI2ME", "Dragen/ EPI2ME", "MEGAN-LR (nucleotide)", "MEGAN-LR (protein)", "Kaiju",
                       "metaMix (nucleotide)", "metaMix Fast (nucleotide)", "metaMix (protein)",  "metaMix Fast (protein)", 
                       "CZ ID* (nucleotide)", "CZ ID* (protein)", "One Codex*")


#### File inputs ####
controls <- read.csv(controls)

community_taxid_mapping <- read.csv(community_taxid_mapping)

genome_lengths <- read.csv(genome_lengths) %>% 
  dplyr::mutate(species = factor(species, msa_2008_levels, msa_2008_labels)) %>% 
  dplyr::rename(name = species)

#### BPM plot ####

bowtie <- read.delim(bowtie_input, sep = "\t", header = TRUE, row.names = NULL)

minimap <- read.delim(minimap_input, sep = "\t", header = TRUE, row.names = NULL)

alignments <- rbind(bowtie, minimap) %>% 
  dplyr::filter(run %in% c(run_levels) & grepl("_sub", sample)) %>%
  dplyr::left_join(controls, by = c("sample", "run")) %>% 
  dplyr::filter(dilution %in% dilutions) %>%
  dplyr::filter(!is.na(control)) %>% 
  dplyr::mutate(run = factor(run, run_levels, run_labels),
         species = factor(species, msa_2008_levels, msa_2008_labels)) %>% 
  dplyr::filter((species %in% c("Human\nbetaherpesvirus 5", "Human\nmastadenovirus F") & dnarna == "DNA" & run != "Twist Comprehensive Viral Panel") |
           (species %in% c("Human\northopneumovirus", "Influenza B\nvirus", "Mammalian\northoreovirus",  "Zika\nvirus") 
            & dnarna == "RNA" & run != "Twist Viral Research Panel (VRP)") |
           (!(species %in% c("Lambda\nphage", "MS2\nphage")) & dnarna == "DNARNA" & run == "Twist Viral Research Panel (VRP)") &
           !is.na(dilution)) %>% 
  dplyr::group_by(run, species, repeat_group, dilution) %>% 
  dplyr::summarise(mean_bp = mean(bp), max_bp = max(bp), min_bp = min(bp)) %>%  
  dplyr::mutate(log10_bp = ifelse(mean_bp == 0, 0, log10(mean_bp)),
         log10_max_bp = ifelse(max_bp == 0, 0, log10(max_bp)),
         log10_min_bp = ifelse(min_bp == 0, 0, log10(min_bp))) %>% 
  dplyr::left_join(genome_lengths) %>% 
  dplyr::mutate(bp_norm = mean_bp*10^4 / genome_length) %>% 
  dplyr::mutate(log10_bp_norm = ifelse(bp_norm == 0, 0, log10(bp_norm)))
  
format_labels <- function(l) {
  
  label <- 6 * 10^l
  
  parse(text=label)
  
}

# Alignment plot
alignments %>% 
  ggplot2::ggplot(aes(x = log10(dilution), y = log10_bp, color = run)) +
  facet_wrap(vars(species)) +
  geom_smooth() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "top") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_x_continuous(breaks = c(1,2,3,4), labels = format_labels) +
  scale_color_manual(values = run_colours) +
  geom_errorbar(aes(ymin = log10_min_bp, max = log10_max_bp), width = 0.1) +
  xlab("Genome copies per ml (gc/ml)") +
  ylab("log10(Bases aligned to genome)")

ggsave(paste0(outdir, "/bp.png"), height = 3.5, width = 6)

# Normalised bases plot
alignments %>% 
  ggplot2::ggplot(aes(x = log10(dilution), y = log10_bp_norm, color = run)) +
  facet_wrap(vars(species)) +
  geom_smooth() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position = "top",
        axis.text.x = element_text(size = 7.5)) +
  scale_x_continuous(breaks = c(1,2,3,4), labels = format_labels) +
  scale_color_manual(values = run_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  xlab("Genome copies per ml (gc/ml)") +
  ylab("Normalised bases")

ggsave(paste0(outdir, "/normalised_bases.png"), height = 3.5, width = 6)

#### All tools sensitivity plot - supplementary ####

empty <- expand.grid(tool = unique(tool_labels_full), species = unique(msa_2008_labels)) %>% 
  merge(controls) %>% 
  filter(run %in% c(run_levels) & grepl("_sub", sample)) %>% 
  mutate(sample_rpm = 0)

report_mc <- read.csv(paste0(cluster_dir, "mock_community_report.csv"))

all_tools <- report_mc %>% 
  filter(run %in% c(run_levels) & grepl("_sub", sample) & tool %in% tool_levels_full) %>% 
  mutate(tool = factor(tool, tool_levels_full, tool_labels_full),
         species = factor(species, msa_2008_levels, msa_2008_labels)) %>%
  select(sample, run, species, tool, control, dnarna, dnarna_pair, repeat_group, repeat., dilution, sample_rpm) %>% 
  rbind(empty) %>% 
  distinct(sample, run, species, tool, control, dnarna, dnarna_pair, repeat_group, repeat., dilution, .keep_all = TRUE) %>%
  group_by(run, species, tool, dnarna, repeat_group, dilution) %>% 
  summarise(mean_sample_rpm = mean(sample_rpm)) %>% 
  mutate(dilution = as.numeric(dilution)) %>% 
  mutate(log10_sample_rpm = ifelse(mean_sample_rpm != 0, log10(mean_sample_rpm*10^6), 0)) %>% 
  filter(!(run %in% c("illumina_250923", "twist_101123") & grepl("megan", tool, ignore.case = TRUE)) & dilution %in% dilutions) %>% 
  mutate(run = factor(run, run_levels, run_labels),
         dilution = as.numeric(dilution)) %>% 
  filter((species %in% c("Human\nbetaherpesvirus 5", "Human\nmastadenovirus F") & dnarna == "DNA" & run != "Twist Viral Research Panel (VRP)") |
           (species %in% c("Human\northopneumovirus", "Influenza B\nvirus", "Mammalian\northoreovirus",  "Zika\nvirus") 
            & dnarna == "RNA" & run != "Twist Viral Research Panel (VRP)") |
           (!(species %in% c("Lambda\nphage", "MS2\nphage")) & dnarna == "DNARNA" & run == "Twist Viral Research Panel (VRP)") &
           !is.na(dilution))

all_tools %>% 
  ggplot2::ggplot(mapping = aes(x = log10(dilution), y = log10_sample_rpm, color = run)) +
  geom_smooth(se = FALSE) +
  facet_grid(rows = vars(species), cols = vars(tool)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_x_continuous(breaks = c(1,2,3,4), labels = format_labels) +
  scale_color_manual(values = run_colours) +
  ylab("log10(reads per million x 10^6)") +
  xlab("Genome copies per ml (gc/ml)")

ggsave(paste0(outdir, "/rpm_alltools.png"), height = 8, width = 13)

