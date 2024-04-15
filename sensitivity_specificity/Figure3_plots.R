library(tidyverse)

thresholds_list <- c(5)

runs <- c("illumina_250923", "nanopore_270923", "twist_101123")

dilutions <- c(10, 100, 1000, 10000)

#### File inputs
controls <- read.csv(paste0(local_dir, "corresponding_controls.csv"))

community_taxid_mapping <- read.csv(paste0(local_dir, "community_taxid_mapping.csv"))

species_counts_input <- paste0(cluster_dir, "species_counts.csv")

false_positives_input <- paste0(cluster_dir, "false_positives_bytype.csv")

#### Levels and labels ####

format_labels <- function(l) {
  
  label <- 6 * 10^l
  
  parse(text=label)
  
}

dnarna_levels <- c("DNA", "RNA", "DNARNA")

dnarna_labels <- c("DNA", "RNA", "DNA+\nRNA")

dilution_levels <- c("10", "100", "1000", "10000")

dilution_labels <- c("60", "600", "6000", "60000")


run_levels <- c("illumina_250923", "nanopore_270923", "twist_101123")

run_labels <- c("Untargeted\nIllumina", "Untargeted\nONT", "Twist Viral\nResearch\nPanel (VRP)")



tool_levels_full <- c("kraken", "bracken", "dragen", "epi2me_bracken", "megan_nucleotide", "megan_protein", "kaiju_refseq-2023-06-08-nucleotide-v2",
                 "metamix_nucleotide", "metamix_step2_nucleotide", "metamix_protein", "metamix_step2_protein", "czid_nt_count", "czid_nr_count", 
                 "onecodex", "onecodex_twist_report")

tool_labels_full <- c("Kraken2", "Bracken", "Dragen/\nEPI2ME", "Dragen/\nEPI2ME", "MEGAN-LR\n(nucleotide)", "MEGAN-LR\n(protein)", "Kaiju",
                 "metaMix\n(nucleotide)", "metaMix Fast\n(nucleotide)", "metaMix \n(protein)",  "metaMix Fast\n(protein)", 
                 "CZ ID*\n(nucleotide)", "CZ ID*\n(protein)", "One Codex*", "One Codex*\n(Twist report)")

tool_levels_main <- c("kraken", "bracken", "dragen", "epi2me_bracken", "megan_nucleotide",
                      "metamix_nucleotide", "metamix_step2_nucleotide", "czid_nt_count", "onecodex", "onecodex_twist_report")

tool_labels_main <- c("Kraken2", "Bracken", "Dragen/\nEPI2ME", "Dragen/\nEPI2ME", "MEGAN-LR",
                 "metaMix", "metaMix Fast", "CZ ID*", "One Codex*", "One Codex*\n(Twist report)")


threshold_levels <- c("0", "1", "2", "5", "9", "10", "11", "12", "13")

threshold_labels <- c("No thresholds", "RPM ratio >= 1", "Set G", "Viruses: PMR > 0.01%, RPMR > 5 or no reads in control\nOther: RPMR > 10, PMR > 1%",
                      "Set A", "Set D", "Set C", "Set E", "Set F")


#### Recall ####

# Empty grid so averages are calculated correctly
empty <- expand.grid(tool = unique(tool_levels_full), thres_ref_no = c(0, thresholds_list)) %>% 
  merge(controls) %>% 
  dplyr::filter(grepl("_sub", sample) & run %in% runs) %>% 
  dplyr::mutate(recall = 0, precision = 0) %>% 
  dplyr::select(run, tool, dnarna_pair, repeat., dilution, thres_ref_no, recall, precision) %>% 
  dplyr::rename(sample = dnarna_pair) %>% 
  dplyr::filter(!(run == "nanopore_270923" & tool == "dragen") &
          !(run %in% c("illumina_250923", "twist_101123") & tool == "epi2me_bracken"))

# This assumes no false positives since we didn't find any with twist onecodex report
onecodex <- read.csv("paste0(local_dir, onecodex_twist_report.csv")) %>% 
  dplyr::select(-c(other_detected, other_indet, n_other_indet)) %>% 
  tidyr::pivot_longer(!c(sample, run, dilution), names_to = "species", values_to = "identified") %>% 
  dplyr::filter(identified == "Y") %>% 
  dplyr::elect(-identified) %>% 
  dplyr::mutate(count = 1,
         dilution = dilution / 6) %>% 
  dplyr::left_join(controls, by = c("sample", "run", "dilution")) %>% 
  dplyr::select(-sample) %>% 
  dplyr::rename(sample = dnarna_pair) %>% 
  dplyr::group_by(run, dilution, sample, repeat.) %>% 
  dplyr::summarise(n_species = sum(count)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(recall = n_species / 6,
         precision = 1,
         tool = "onecodex_twist_report",
         thres_ref_no = 0) %>% 
  dplyr::select(-n_species)

species_counts_plot <- read.csv(species_counts_input) %>% 
  dplyr::rename(sample = dnarna_pair) %>%
  dplyr::filter(tool %in% tool_levels_full & thres_ref_no %in% c(0, thresholds_list)) %>% 
  dplyr::filter(grepl("_sub", sample) & run %in% runs) %>% 
  dplyr::mutate(dilution = as.numeric(dilution),
          thres_ref_no = as.character(thres_ref_no)) %>% 
  dplyr::select(run, tool, sample,  repeat., dilution, recall, precision, thres_ref_no) %>%
  rbind(onecodex) %>% 
  rbind(empty) %>%
  dplyr::arrange(desc(recall)) %>% 
  dplyr::distinct(run, tool, sample, repeat., dilution, thres_ref_no, .keep_all = TRUE) %>%
  dplyr::filter(dilution %in% dilutions) %>% 
  dplyr::group_by(run, tool, dilution, thres_ref_no) %>% 
  dplyr::summarise(mean_recall = mean(recall), mean_precision = mean(precision))

# Plots coomparing recall at different thresholds '
recall_plot <- function(species_counts_plot, tool_levels, tool_labels, thresholds_list, outdir, name, height, width) {
  
  for (thres_ref in thresholds_list) {
  
    plot <- species_counts_plot %>% 
      dplyr::filter(thres_ref_no %in% c(thres_ref, 0)) %>% 
      dplyr::filter(tool %in% tool_levels) %>% 
      dplyr::mutate(tool = factor(tool, tool_levels, tool_labels),
             run = factor(run, run_levels, run_labels),
             dilution = factor(as.character(dilution), dilution_levels, dilution_labels),
             thres_ref_no = factor(thres_ref_no, threshold_levels, threshold_labels)) %>% 
      ggplot2::ggplot(aes(x = dilution, y = mean_recall, fill = thres_ref_no)) +
      geom_col(position = position_dodge2(preserve = "single"), width = 0.6) +
      facet_grid(cols = vars(tool), rows = vars(run)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            legend.position = "top",
            legend.title = element_blank(),
            strip.text.y = element_text(angle = 0),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_y_continuous(expand = expansion(mult = c(0, .1))) +
      ylab("Sensitivity") +
      xlab("Genome copies per ml (gc/ml)")+
      ggtitle(paste0("Thresholds: ", thres_ref))
    
    recall_output <- paste0(outdir, "/recall_", name, "_thres", thres_ref, ".png")
  
    ggsave(recall_output, plot, height = height, width = width)
    #return(plot)
    
  }

}

recall_plot(species_counts_plot, tool_levels_full, tool_labels_full, thresholds_list, outdir, name = "full", height = 4.5, width = 13)
recall_plot(species_counts_plot, tool_levels_main, tool_labels_main, thresholds_list, outdir, name = "main", height = 4.5, width = 10)

recall_plot(filter(species_counts_plot, tool != "onecodex_twist_report"), tool_levels_main, tool_labels_main, 
                     thresholds_list, outdir, name = "SI", height = 4.5, width = 10)


#### Specificity ####

onecodex_specificity <- read.csv(paste0(local_dir, "onecodex_twist_report_fp.csv")) %>% 
  dplyr::left_join(controls) %>% 
  dplyr::mutate(db = NA, type = "Virus", tool = "onecodex_twist_report", thres_ref_no = 0) %>% 
  dplyr::select(-sample) %>% 
  dplyr::rename(sample = dnarna_pair) %>% 
  dplyr::select(sample, run, db, repeat., dilution, type, tool, thres_ref_no, false_positives)

false_positives_bytype <- read.csv(false_positives_input) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>%
  dplyr::rename(sample = dnarna_pair) %>% 
  dplyr::filter(tool %in% tool_levels_full & grepl("_sub", sample) & run %in% runs & dilution %in% dilutions) %>% 
  rbind(onecodex_specificity) %>% 
  dplyr::group_by(run, tool, type, dilution, thres_ref_no) %>% 
  dplyr::summarise(mean_false_positives = mean(false_positives)) %>% 
  dplyr::mutate(dilution = as.numeric(dilution))
 
type_colours <- setNames(c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
                         c("Bacteria", "Fungi", "Other eukaryote", "Virus"))

false_positives_plot <- function(false_positives_bytype, tool_levels, tool_labels, thresholds_list, outdir, name, height, width) {

  for (thres_ref in c(0, thresholds_list)) {
    
    plot <- false_positives_bytype %>%
      dplyr::filter(thres_ref_no == thres_ref & tool %in% tool_levels) %>% 
      dplyr::mutate(tool = factor(tool, tool_levels, tool_labels),
             run = factor(run, run_levels, run_labels),
             dilution = factor(as.character(dilution), dilution_levels, dilution_labels),
             thres_ref_no = factor(thres_ref_no, threshold_levels, threshold_labels)) %>%
      ggplot(aes(x = dilution, y = mean_false_positives, fill = type)) +
      geom_bar(stat = "identity") +
      facet_grid(cols = vars(tool), rows = vars(run)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            legend.position = "top",
            legend.title = element_blank(),
            strip.text.y = element_text(angle = 0),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_y_continuous(expand = expansion(mult = c(0, .1))) +
      scale_fill_manual(values = type_colours) +
      ylab("Number of false positive species") +
      xlab("Genome copies per ml (gc/ml)") +
      ggtitle(paste0("Thresholds: ", thres_ref))
    
    false_positives_plot_output <- paste0(outdir, "/false_positives_", name, "_thres", thres_ref, ".png")
    
    ggsave(false_positives_plot_output, plot, height = height, width = width)
    #return(plot)
    
  }

}

false_positives_plot(false_positives_bytype, tool_levels_full, tool_labels_full, thresholds_list, outdir, name = "full", height = 4.5, width = 13)
false_positives_plot(false_positives_bytype, tool_levels_main, tool_labels_main, thresholds_list, outdir, name = "main", height = 4.5, width = 10)

false_positives_plot(filter(false_positives_bytype, tool != "onecodex_twist_report"), tool_levels_main, tool_labels_main, thresholds_list, outdir, name = "SI", height = 4.5, width = 10)


# Viruses only  
false_positives_plot_viruses <- function(false_positives_bytype, tool_levels, tool_labels, thresholds_list, outdir, name, height, width, bar_width) {

  for (thres_ref in thresholds_list) { 
  
    plot <- false_positives_bytype %>% 
      dplyr::filter(thres_ref_no %in% c(0, thres_ref) & type == "Virus" & tool %in% tool_levels) %>% 
      dplyr::mutate(tool = factor(tool, tool_levels, tool_labels),
             run = factor(run, run_levels, run_labels),
             dilution = factor(as.character(dilution), dilution_levels, dilution_labels),
             thres_ref_no = factor(thres_ref_no, threshold_levels, threshold_labels)) %>% 
      ggplot2::ggplot(aes(x = dilution, y = mean_false_positives, fill = thres_ref_no)) +
      geom_col(position = position_dodge2(preserve = "single"), width = bar_width) +
      facet_grid(cols = vars(tool), rows = vars(run)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            legend.position = "top",
            legend.title = element_blank(),
            strip.text.y = element_text(angle = 0),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_y_continuous(expand = expansion(mult = c(0, .1))) +
      ylab("Number of false positive viral species") +
      xlab("Genome copies per ml (gc/ml)") +
      ggtitle(paste0("Thresholds: ", thres_ref))
    
    false_positives_plot_output <- paste0(outdir, "/false_positive_viruses_", name, "_thres", thres_ref, ".png")
    
    ggsave(false_positives_plot_output, height = height, width = width)
    #return(plot)
    
  }
  
}

false_positives_plot_viruses(false_positives_bytype, tool_levels_full, tool_labels_full, thresholds_list, outdir, name = "full", height = 4.5, width = 13, bar_width = 0.6)
false_positives_plot_viruses(false_positives_bytype, tool_levels_main, tool_labels_main, thresholds_list, outdir, name = "main", height = 4.5, width = 10, bar_width = 0.6)

false_positives_plot_viruses(filter(false_positives_bytype, tool != "onecodex_twist_report"), tool_levels_main, tool_labels_main, 0, outdir, 
                             name = "SI", height = 4.5, width = 10, bar_width = 0.9)
false_positives_plot_viruses(filter(false_positives_bytype, tool != "onecodex_twist_report"), tool_levels_main, tool_labels_main, thresholds_list, outdir, 
                             name = "SI", height = 4.5, width = 10, bar_width = 0.6)


# Further analysis of false positive species
msa_2008_levels <- c("Human betaherpesvirus 5", "Human herpesvirus 5", "Cytomegalovirus humanbeta5",
                     "Human mastadenovirus F", 
                     "Human orthopneumovirus", "Orthopneumovirus hominis", "Human respiratory syncytial virus",
                     "Influenza B virus", "Betainfluenzavirus influenzae",
                     "Mammalian orthoreovirus", "Reovirus 3",
                     "Zika virus", "Orthoflavivirus zikaense",
                     "Lambda phage", "Escherichia virus Lambda", "Lambdavirus lambda",
                     "MS2 phage", "Escherichia virus MS2", "Emesvirus zinderi")



positive_species <- read.csv(paste0(cluster_dir, "positive_species.csv"))

false_positive_viruses <- positive_species %>% 
  dplyr::filter(type == "Virus" & thres_ref_no == 5 & tool %in% tool_levels_main & !(species %in% msa_2008_levels) & 
           grepl("_sub", sample) & run %in% run_levels) %>% 
  dplyr::select(run, dilution, tool, species, species_taxid) %>% 
  dplyr::distinct() %>% 
  dplyr::arrange(run, dilution, tool, species)

fp_virus_list <- false_positive_viruses %>% 
  dplyr::select(species, species_taxid) %>% 
  dplyr::arrange(species) %>% 
  dplyr::distinct()

write.csv(fp_virus_list, paste0(cluster_dir, "false_positive_virus_types2.csv"), quote = FALSE, row.names = FALSE)

fp_virus_types <- read.csv(paste0(cluster_dir, "false_positive_virus_types.csv"))

virus_class <- function(infects, virus_type) {
  
  if (virus_type == "Related") {
    
    return("related")
    
  } else if (virus_type == "Anelloviridae") {
    
    return("Anelloviridae")
    
  } else if (infects %in% c("humans", "mammals", "birds")) {
    
    return("humans_mammals_birds_other")
    
  } else {
    
    return(infects)
    
  }
  
}

virus_class <- Vectorize(virus_class)

false_positive_viruses2 <- false_positive_viruses %>% 
  dplyr::left_join(fp_virus_types, by = c("species", "species_taxid")) %>% 
  dplyr::mutate(virus_class = virus_class(infects, virus_type), count = 1) %>% 
  dplyr::group_by(run, dilution, tool, virus_class) %>% 
  dplyr::summarise(n_species = sum(count))

virus_class_levels <- c("related", "Anelloviridae", "humans_mammals_birds_other",
                        "fungi", "invertebrates", "protists", "fishes", "plants", 
                        "bacteria", "archaea")

virus_class_labels <- c("Related to virus in mock community", "Anelloviridae", "Other human/mammalian viruses",
                        "Infects other eukaryote", "Infects other eukaryote", "Infects other eukaryote", "Infects other eukaryote", "Infects other eukaryote",
                        "Phage", "Phage")
  
  
false_positive_viruses2 %>% 
  dplyr::mutate(tool = factor(tool, tool_levels_main, tool_labels_main),
         run = factor(run, run_levels, run_labels),
         virus_class = factor(virus_class, virus_class_levels, virus_class_labels)) %>%
  ggplot2::ggplot(aes(x = log10(dilution), y = n_species, fill = virus_class)) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(tool), rows = vars(run)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_continuous(breaks = c(1,2,3,4), labels = format_labels) +
  ylab("Number of distinct false positive viral species") +
  xlab("Genome copies per ml (gc/ml)")

ggsave(paste0(outdir, "/false_positive_virus_class_thres5.png"), height = 4.5, width = 10)

