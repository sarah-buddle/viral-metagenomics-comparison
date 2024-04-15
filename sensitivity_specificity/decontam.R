library(tidyverse)
library(phyloseq)
library(decontam)

taxonomy_csv <- ""
report_input <- ""
controls <- ""

runs <- c("illumina_250923", "nanopore_270923", "twist_101123")

tool_levels_main <- c("kraken", "bracken", "dragen", "epi2me_bracken", "megan_nucleotide",
                      "metamix_nucleotide", "metamix_step2_nucleotide", "czid_nt_count", "onecodex")

report <- read.csv(report_input) %>% 
  dplyr::mutate(db = replace_na(db, "NA")) %>% 
  dplyr::filter(!is.na(species_taxid)) %>% 
  dplyr::filter(run %in% runs & tool %in% tool_levels_main & grepl("_sub", sample) & !grepl("sub2|sub5", sample)) %>% 
  dplyr::filter(!(run == "twist_101123" & tool == "onecodex"))

# Make OTU matrix
otumat <- report %>% 
  dplyr::mutate(name = paste0(sample, "_", run, "_", tool)) %>% 
  dplyr::select(name, species_taxid, reads) %>% 
  tidyr::pivot_wider(names_from = "name", values_from = "reads") %>% 
  tibble::column_to_rownames(var = "species_taxid") %>% 
  dplyr::mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  as.matrix(.)

otu <- phyloseq::otu_table(otumat, taxa_are_rows = TRUE)

# Make taxonomy matrix
taxonomy <- read.csv(taxonomy_csv) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(taxid = as.numeric(taxid), species_taxid = as.numeric(species_taxid))

taxmat <- taxonomy %>% 
  dplyr::select(species_taxid, type, species) %>% 
  dplyr::filter(!is.na(species_taxid)) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(kingdom = type) %>% 
  tibble::column_to_rownames(var = "species_taxid") %>% 
  as.matrix(.)

tax <- phyloseq::tax_table(taxmat)

# Sample metadata
tools <- expand_grid(tool = tool_levels_main, run = runs)

sampledf <- read.csv(controls) %>% 
  dplyr::filter(run %in% runs & grepl("_sub", sample) & !grepl("sub2|sub5", sample)) %>% 
  dplyr::full_join(tools, by = "run") %>% 
  dplyr::mutate(name = paste0(sample, "_", run, "_", tool)) %>% 
  dplyr::select(-c(sample, run, tool)) %>% 
  tibble::column_to_rownames(var = "name") %>% 
  dplyr::mutate(conc_ngperul = 40) %>% 
  dplyr::mutate(control = gsub("control", "Control Sample", control)) %>% 
  dplyr::mutate(control = gsub("sample", "True Sample", control)) %>% 
  dplyr::rename(Sample_or_Control = control)

samples <- phyloseq::sample_data(sampledf)

ps <- phyloseq::phyloseq(otu, samples, tax)

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control Sample"

community_taxid_mapping <- read.csv(community_taxid_mapping)

phage_taxids <- c(329852, 2681611, 10710)

roc <- data.frame()

# Make roc curve by varying thresholds
for (threshold in c(seq(0,1,0.01), c(0.999, 0.9999, 0.99999, 0.999999, 0.9999999, 0.99999999, 0.999999999, 0.9999999999, 0.99999999999))) {

  contamdf.freq <- isContaminant(ps, method="prevalence", neg = "is.neg", conc="conc_ngperul", threshold = threshold)
  
  decontam_df <- contamdf.freq %>% 
    dplyr::select(contaminant) %>% 
    tibble::rownames_to_column(var = "species_taxid") %>% 
    dplyr::mutate(species_taxid = as.numeric(species_taxid))
  
  report_decontam <- report %>% 
    dplyr::filter(!(species_taxid %in% phage_taxids)) %>% 
    dplyr::left_join(decontam_df) %>% 
    dplyr::filter(type == "Virus")
  
  species <- unique(report_decontam$species_taxid)
  
  report_decontam_positive <- report_decontam %>% 
    dplyr::filter(contaminant == "FALSE")
  
  positive_species <- unique(report_decontam_positive$species_taxid)
  
  expected_taxids <- community_taxid_mapping %>% 
    dplyr::filter(community == "msa_2008" & positive_control == "FALSE")
  
  expected_true_positives <- unique(expected_taxids$taxid)
  
  true_positives <- as.numeric(sum(expected_true_positives %in% positive_species))
  
  false_positives <- as.numeric(sum(!(positive_species %in% expected_true_positives)))
  
  false_negatives <- length(expected_true_positives) - true_positives
  
  true_negatives <- length(species) - true_positives - false_positives - false_negatives
  
  new_row <- c(threshold, true_positives, false_positives, false_negatives, true_negatives)
  
  roc <- rbind(roc, new_row)
}

colnames(roc) <- c("threshold", "true_positives", "false_positives", "false_negatives", "true_negatives")

roc_ss <- roc %>% 
  mutate(sensitivity = true_positives / (true_positives + false_negatives),
         specificity = true_negatives / (true_negatives + false_positives))

write.csv(roc_ss, output, quote = FALSE, row.names = FALSE)
  
