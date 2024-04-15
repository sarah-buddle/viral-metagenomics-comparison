library(tidyverse)

options(scipen=999)

# Mappings  
community_taxid_mapping <- read.csv(paste0(local_dir, "community_taxid_mapping.csv"))

run_community_mapping <- read.csv(paste0(local_dir, "run_community_mapping.csv"))

mappings <- full_join(community_taxid_mapping, run_community_mapping, by = "community") %>% 
  dplyr::rename(species = name)

controls <- read.csv(paste0(local_dir, "corresponding_controls_withcontrols.csv"))

# Read in CZID 
czid_raw_nanopore <- read.csv(paste0(cluster_dir, "czid/all_czid_nanopore.csv"))

czid_raw_illumina <- read.csv(paste0(cluster_dir, "czid/all_czid_illumina.csv"))

czid_common_cols <- intersect(colnames(czid_raw_nanopore), colnames(czid_raw_illumina))

czid_raw <- full_join(czid_raw_illumina, czid_raw_nanopore, by = czid_common_cols)

czid <- czid_raw %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>%
  dplyr::rename(taxid = tax_id) %>% 
  dplyr::select(c(sample, run, name, taxid, nt_count, nr_count)) %>% 
  tidyr::pivot_longer(!c(sample, run, name, taxid), names_to = "db", values_to = "reads") %>% 
  dplyr::mutate(tool = "czid") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = TRUE)

# Preprocess kraken

kraken_raw <- read.delim(paste0(cluster_dir, "kraken2/all_kraken.txt"), sep = "\t", header = TRUE, row.names = NULL)

kraken <- kraken_raw %>%
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(tool = "kraken") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = TRUE)


# Preprocess_bracken

bracken_raw <- read.delim(paste0(cluster_dir, "bracken/all_bracken.txt"), sep = "\t", header = TRUE, row.names = NULL)

bracken <- bracken_raw %>% 
  dplyr::select(-c("taxonomy_lvl", "kraken_assigned_reads", "added_reads", "fraction_total_reads")) %>% 
  dplyr::rename(taxid = taxonomy_id, reads = new_est_reads) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>%
  dplyr::mutate(tool = "bracken") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = TRUE)


# Preprocess dragen

dragen_input <- paste0(cluster_dir, "dragen/all_dragen.txt")

dragen_raw <- read.delim(dragen_input, sep = "\t", header = TRUE, row.names = NULL)

dragen <- dragen_raw %>% 
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(tool = "dragen") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = TRUE)

# Preprocess epi2me bracken

epi2me_bracken_raw <- read.delim(paste0(cluster_dir, "epi2me_kraken/all_epi2me_bracken.txt"), 
                                 sep = "\t", header = TRUE, row.names = NULL)

epi2me_bracken <- epi2me_bracken_raw %>% 
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(tool = "epi2me_bracken") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = TRUE)

# Preprocess metamix

metamix_raw <- read.delim(paste0(cluster_dir, "metamix/all_metamix.txt"), sep = "\t", header = TRUE, row.names = NULL)

metamix <- metamix_raw %>% 
  dplyr::filter(name != "unknown" & log10BF > 0) %>% 
  dplyr::select(c(sample, run, db, name, taxid, reads)) %>%
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(tool = "metamix") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = TRUE)

# Preprocess metaMix step 2

metamix_step2_raw <- read.delim( paste0(cluster_dir, "metamix/all_metamix_step2.txt"), sep = "\t", header = TRUE, row.names = NULL)

metamix_step2 <- metamix_step2_raw %>% 
  dplyr::select(-c(rowNumber, samplingWeight)) %>% 
  dplyr::rename(name = scientName, taxid = taxonID, reads = countReads) %>%
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(tool = "metamix_step2") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = TRUE)


# Preprocess Megan

megan_raw <- read.delim(paste0(cluster_dir, "megan/all_megan.txt"), sep = "\t", header = TRUE, row.names = NULL)

megan <- megan_raw %>% 
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(tool = "megan") %>% 
  dplyr::mutate(tool = ifelse(grepl("nucleotide", db), "megan_nucleotide", "megan_protein")) %>% 
  dplyr::select(-db)

# Preprocess kaiju

kaiju_raw <- read.delim(paste0(cluster_dir, "kaiju/all_kaiju.txt"), sep = "\t", header = TRUE, row.names = NULL)

kaiju <- kaiju_raw %>% 
  dplyr::select(-c(file, percent)) %>% 
  dplyr::rename(name = taxon_name) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(tool = "kaiju") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = TRUE)

# Preprocess bowtie

bowtie_raw <- read.delim(paste0(cluster_dir, "bowtie/all_bowtie.txt"),
                         sep = "\t", header = TRUE, row.names = NULL)

bowtie <- bowtie_raw %>% 
  dplyr::mutate(species = gsub(" ", "_", species)) %>%
  dplyr::mutate(species = gsub("_reference", "", species)) %>%
  dplyr::left_join(mappings, by = "species") %>% 
  dplyr::mutate(species = gsub("_", " ", species)) %>%
  dplyr::select(-c(taxid_czid, bp, community, run.y, positive_control, platform)) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(run = run.x, name = species)


bowtie_nodedup_raw <- read.delim(paste0(cluster_dir, "bowtie/all_bowtie_nodedup.txt"),
                                 sep = "\t", header = TRUE, row.names = NULL)

bowtie_nodedup <- bowtie_nodedup_raw %>% 
  dplyr::mutate(species = gsub(" ", "_", species)) %>%
  dplyr::mutate(species = gsub("_reference", "", species)) %>%
  dplyr::left_join(mappings, by = "species") %>% 
  dplyr::mutate(species = gsub("_", " ", species)) %>%
  dplyr::select(-c(taxid_czid, bp, community, run.y, positive_control, platform)) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(run = run.x, name = species) %>% 
  dplyr::mutate(tool = "bowtie_nodedup")

# Preprocess minimap

minimap_raw <- read.delim(paste0(cluster_dir, "minimap/all_minimap.txt"),
                          sep = "\t", header = TRUE, row.names = NULL)

minimap <- minimap_raw %>% 
  dplyr::mutate(species = gsub(" ", "_", species)) %>%
  dplyr::mutate(species = gsub("_reference", "", species)) %>%
  dplyr::left_join(mappings, by = "species") %>% 
  dplyr::mutate(species = gsub("_", " ", species)) %>%
  dplyr::select(-c(taxid_czid, bp, community, run.y, positive_control, platform)) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(run = run.x, name = species)

minimap_nodedup_raw <- read.delim(paste0(cluster_dir, "minimap/all_minimap_nodedup.txt"), 
                                   sep = "\t", header = TRUE, row.names = NULL)

minimap_nodedup <- minimap_nodedup_raw %>% 
  dplyr::mutate(species = gsub(" ", "_", species)) %>%
  dplyr::mutate(species = gsub("_reference", "", species)) %>%
  dplyr::left_join(mappings, by = "species") %>% 
  dplyr::mutate(species = gsub("_", " ", species)) %>%
  dplyr::select(-c(taxid_czid, bp, community, run.y, positive_control, platform)) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(run = run.x, name = species) %>% 
  dplyr::mutate(tool = "minimap_nodedup")

# Onecodex

onecodex <- read.csv(paste0(cluster_dir, "onecodex/onecodex.csv")) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>%
  dplyr::select(Sample.Name, Tax.Name, Tax.ID, Reads) %>%
  dplyr::rename(samplerun = Sample.Name, name = Tax.Name, taxid = Tax.ID, reads = Reads) %>% 
  dplyr::mutate(samplerun = gsub(".fastq.gz", "", samplerun)) %>% 
  dplyr::mutate(samplerun = gsub(".fq.gz", "", samplerun)) %>% 
  dplyr::mutate(run = str_extract(samplerun, "([^_\n]+)_[^_]*$"),
         sample = gsub("_twist_101123|_illumina_250923|_nanopore_270923", "", samplerun)) %>% # change this part if possible
  dplyr::select(-samplerun) %>% 
  dplyr::mutate(tool = "onecodex")

# Levels and labels

dilution_levels <- c("10", "100", "1000", "10000")

dilution_labels <- c("60", "600", "6000", "60000")


run_levels <- c("illumina_250923", "nanopore_270923", "twist_101123")

run_labels <- c("Untargeted Illumina", "Untargeted ONT", "Twist VRP")



tool_levels_full <- c("bowtie", "bowtie_nodedup", "minimap",  "minimap_nodedup",
                      "kraken_refseq-2023-06-08-nucleotide-v2", "bracken_refseq-2023-06-08-nucleotide-v2",
                      "dragen_refseq-2023-06-08-nucleotide-v2", "epi2me_bracken_refseq-2023-06-08-nucleotide-v2",
                      "megan_nucleotide", "megan_protein", "kaiju_refseq-2023-06-08-nucleotide-v2",
                      "metamix_nucleotide", "metamix_step2_nucleotide", "metamix_protein", "metamix_step2_protein",
                      "czid_nt_count", "czid_nr_count", 
                      "onecodex", "onecodex_twist_report")

tool_labels_full <- c("Bowtie2 (deduplicated)", "Bowtie2 (no deduplication)", "Minimap2 (deduplicated)", "Minimap2 (no deduplication)",
                      "Kraken2", "Bracken", 
                      "Dragen", "EPI2ME",
                      "MEGAN-LR (nucleotide)", "MEGAN-LR (protein)", "Kaiju",
                      "metaMix (nucleotide)", "metaMix Fast (nucleotide)", "metaMix (protein)",  "metaMix Fast (protein)", 
                      "CZ ID* (nucleotide)", "CZ ID* (protein)",
                      "One Codex*", "One Codex* (Twist report)")




# Combine
report <- rbind(czid, dragen, epi2me_bracken, kraken, bracken, metamix, metamix_step2, bowtie, bowtie_nodedup, minimap, minimap_nodedup, megan, kaiju, onecodex) %>%  
  dplyr::distinct() %>% 
  dplyr::left_join(controls) %>% 
  dplyr::filter(!is.na(name) & !is.na(reads) & reads != 0 & 
                  run %in% run_levels & dilution %in% dilution_levels & grepl("sub", sample) & tool %in% tool_levels_full) %>% 
  dplyr::select(run, dnarna, repeat., dilution, tool, name, taxid, reads) %>% 
  dplyr::mutate(run = factor(run, run_levels, run_labels),
                tool = factor(tool, tool_levels_full, tool_labels_full),
                dilution = factor(dilution, dilution_levels, dilution_labels)) %>% 
  dplyr::arrange(run, dnarna, repeat., desc(dilution), tool, desc(reads)) %>%
  dplyr::rename(platform = run, gcperml = dilution)

write.csv(report, paste0(cluster_dir, "raw_report_alltaxa.csv"), row.names = FALSE, quote = FALSE)


