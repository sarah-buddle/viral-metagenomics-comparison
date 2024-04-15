library(tidyverse)

options(scipen=999)

# Mappings  
community_taxid_mapping <- read.csv(paste0(local_dir, "community_taxid_mapping.csv"))

run_community_mapping <- read.csv(paste0(local_dir, "run_community_mapping.csv"))

mappings <- full_join(community_taxid_mapping, run_community_mapping, by = "community") %>% 
  dplyr::rename(species = name)

controls <- read.csv(paste0(local_dir, "corresponding_controls.csv"))

# Read counts

read_counts <- read.csv(paste0(cluster_dir, "counts.csv")) %>%  
  select(c(sample, run, raw_reads))

# Read in taxonomy
taxonomy <- read.csv(taxonomy_csv) %>% 
  distinct() %>% 
  mutate(taxid = as.numeric(taxid), species_taxid = as.numeric(species_taxid))


# Be sure to remove quotes from file before import
czid_taxonomy <- read.csv(czid_taxonomy_csv, quote = "") %>% 
  mutate(taxid = as.numeric(taxid), species_taxid = as.numeric(species_taxid))

# Read in CZID 
czid_raw_nanopore <- read.csv(paste0(cluster_dir, "czid/all_czid_nanopore.csv"))

czid_raw_illumina <- read.csv(paste0(cluster_dir, "czid/all_czid_illumina.csv"))

czid_common_cols <- intersect(colnames(czid_raw_nanopore), colnames(czid_raw_illumina))

czid_raw <- full_join(czid_raw_illumina, czid_raw_nanopore, by = czid_common_cols)

czid <- czid_raw %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>%
  dplyr::filter(tax_level == 1) %>% 
  dplyr::rename(taxid = tax_id) %>% 
  dplyr::select(c(sample, run, taxid, nt_count, nr_count)) %>% 
  dplyr::left_join(czid_taxonomy, by = "taxid") %>% 
  dplyr::mutate(species = ifelse(taxid == -100, "unclassified", species)) %>% 
  dplyr::select(!c(taxid, name)) %>% 
  dplyr::group_by(sample, run, species, species_taxid, type) %>% 
  dplyr::summarise(nt_count = sum(nt_count, na.rm = TRUE), nr_count = sum(nr_count, na.rm = TRUE)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(!c(sample, run, species, species_taxid, type), names_to = "db", values_to = "reads") %>% 
  dplyr::mutate(tool = "czid") %>% 
  dplyr::filter(reads != 0) %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = FALSE)

# Preprocess kraken

kraken_raw <- read.delim(paste0(cluster_dir, "kraken2/all_kraken.txt"), sep = "\t", header = TRUE, row.names = NULL)

kraken <- kraken_raw %>%
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "kraken")


# Preprocess_bracken

bracken_raw <- read.delim(paste0(cluster_dir, "bracken/all_bracken.txt"), sep = "\t", header = TRUE, row.names = NULL)

bracken <- bracken_raw %>% 
  dplyr::select(-c("taxonomy_lvl", "name", "kraken_assigned_reads", "added_reads", "fraction_total_reads")) %>% 
  dplyr::rename(taxid = taxonomy_id, reads = new_est_reads) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(taxid = as.numeric(taxid), reads = as.numeric(reads)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "bracken")

# Preprocess dragen

dragen_input <- paste0(cluster_dir, "dragen/all_dragen.txt")

dragen_raw <- read.delim(dragen_input, sep = "\t", header = TRUE, row.names = NULL)

dragen <- dragen_raw %>% 
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "dragen")

# Preprocess epi2me kraken

epi2me_kraken_raw <- read.delim(paste0(cluster_dir, "epi2me_kraken/all_epi2me_kraken.txt"), sep = "\t", header = TRUE, row.names = NULL)

epi2me_kraken <- epi2me_kraken_raw %>% 
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "epi2me_kraken")

# Preprocess epi2me bracken

epi2me_bracken_raw <- read.delim(paste0(cluster_dir, "epi2me_kraken/all_epi2me_bracken.txt"), 
                                 sep = "\t", header = TRUE, row.names = NULL)

epi2me_bracken <- epi2me_bracken_raw %>% 
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "epi2me_bracken")

# Preprocess metamix

metamix_raw <- read.delim(paste0(cluster_dir, "metamix/all_metamix.txt"), sep = "\t", header = TRUE, row.names = NULL)

metamix <- metamix_raw %>% 
  dplyr::filter(name != "unknown" & log10BF > 0) %>% 
  dplyr::select(c(sample, run, db, name, taxid, reads)) %>%
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(taxid = as.numeric(taxid)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "metamix") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = FALSE)

# Preprocess metaMix step 2

metamix_step2_raw <- read.delim( paste0(cluster_dir, "metamix/all_metamix_step2.txt"), sep = "\t", header = TRUE, row.names = NULL)

metamix_step2 <- metamix_step2_raw %>% 
  dplyr::select(-c(rowNumber, samplingWeight)) %>% 
  dplyr::rename(name = scientName, taxid = taxonID, reads = countReads) %>%
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::mutate(taxid = as.numeric(taxid)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "metamix_step2") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = FALSE)


# Preprocess Megan

megan_raw <- read.delim(paste0(cluster_dir, "megan/all_megan.txt"), sep = "\t", header = TRUE, row.names = NULL)

megan <- megan_raw %>% 
  dplyr::select(c(sample, run, db, name, taxid, assignedHere)) %>% 
  dplyr::rename(reads = assignedHere) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "megan") %>% 
  dplyr::mutate(tool = ifelse(grepl("nucleotide", db), "megan_nucleotide", "megan_protein"))

# Preprocess kaiju

kaiju_raw <- read.delim(paste0(cluster_dir, "kaiju/all_kaiju.txt"), sep = "\t", header = TRUE, row.names = NULL)

kaiju <- kaiju_raw %>% 
  dplyr::select(-c(file, percent, taxon_name)) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "kaiju") %>% 
  tidyr::unite(tool, c(tool, db), sep = "_", remove = FALSE)

# Preprocess bowtie

bowtie_raw <- read.delim(paste0(cluster_dir, "bowtie/all_bowtie.txt"),
                         sep = "\t", header = TRUE, row.names = NULL)

bowtie <- bowtie_raw %>% 
  dplyr::mutate(db = NA) %>% 
  dplyr::mutate(species = gsub(" ", "_", species)) %>%
  dplyr::mutate(species = gsub("_reference", "", species)) %>%
  dplyr::left_join(mappings, by = "species") %>% 
  dplyr::select(-c(taxid_czid, bp, community, run.y)) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(species = gsub("_", " ", species)) %>%
  dplyr::rename(run = run.x, name = species) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type, tool) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::filter(species != "")


bowtie_nodedup_raw <- read.delim(paste0(cluster_dir, "bowtie/all_bowtie_nodedup.txt"),
                                 sep = "\t", header = TRUE, row.names = NULL)

bowtie_nodedup <- bowtie_nodedup_raw %>% 
  dplyr::mutate(db = NA) %>% 
  dplyr::mutate(species = gsub(" ", "_", species)) %>%
  dplyr::mutate(species = gsub("_reference", "", species)) %>%
  dplyr::left_join(mappings, by = "species") %>% 
  dplyr::select(-c(taxid_czid, bp, community, run.y)) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(species = gsub("_", " ", species)) %>%
  dplyr::rename(run = run.x, name = species) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type, tool) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::filter(species != "") %>% 
  dplyr::mutate(tool = "bowtie_nodedup")

# Preprocess minimap

minimap_raw <- read.delim(paste0(cluster_dir, "minimap/all_minimap.txt"),
                          sep = "\t", header = TRUE, row.names = NULL)

minimap <- minimap_raw %>% 
  dplyr::mutate(db = NA) %>% 
  dplyr::mutate(species = gsub(" ", "_", species)) %>%
  dplyr::mutate(species = gsub("_reference", "", species)) %>%
  dplyr::left_join(mappings, by = "species") %>% 
  dplyr::select(-c(taxid_czid, bp, community, run.y)) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(run = run.x, name = species) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type, tool) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::filter(!is.na(species_taxid))

minimap_nodedup_raw <- read.delim(paste0(cluster_dir, "minimap/all_minimap_nodedup.txt"), 
                                   sep = "\t", header = TRUE, row.names = NULL)

minimap_nodedup <- minimap_nodedup_raw %>% 
  dplyr::mutate(db = NA) %>% 
  dplyr::mutate(species = gsub(" ", "_", species)) %>%
  dplyr::mutate(species = gsub("_reference", "", species)) %>%
  dplyr::left_join(mappings, by = "species") %>% 
  dplyr::select(-c(taxid_czid, bp, community, run.y)) %>% 
  dplyr::distinct() %>% 
  dplyr::rename(run = run.x, name = species) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type, tool) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::filter(!is.na(species_taxid)) %>% 
  dplyr::mutate(tool = "minimap_nodedup")

# Onecodex

onecodex <- read.csv(paste0(cluster_dir, "onecodex/onecodex.csv")) %>% 
  dplyr::mutate(across(where(is.character), str_trim)) %>%
  dplyr::select(Sample.Name, Tax.ID, Reads) %>%
  dplyr::rename(samplerun = Sample.Name, taxid = Tax.ID, reads = Reads) %>% 
  dplyr::mutate(samplerun = gsub(".fastq.gz", "", samplerun)) %>% 
  dplyr::mutate(samplerun = gsub(".fq.gz", "", samplerun)) %>% 
  dplyr::mutate(db = NA,
         run = str_extract(samplerun, "([^_\n]+)_[^_]*$"),
         sample = gsub("_twist_101123|_illumina_250923|_nanopore_270923", "", samplerun)) %>% # change this part if possible
  dplyr::select(-samplerun) %>% 
  dplyr::left_join(taxonomy, by = "taxid") %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type) %>% 
  dplyr::summarise(reads = sum(reads, na.rm = TRUE)) %>% 
  dplyr::mutate(tool = "onecodex")


# Combiens reports 
report <- rbind(czid, dragen, epi2me_bracken, kraken, bracken, metamix, metamix_step2, bowtie, bowtie_nodedup, minimap, minimap_nodedup, megan, kaiju, onecodex) %>%  
  dplyr::distinct() %>% 
  dplyr::filter(is.na(species) | species != "Homo sapiens") %>% 
  dplyr::filter(reads != 0) %>% 
  adplyr::rrange(type)

write.csv(report, paste0(cluster_dir, "raw_report.csv"), row.names = FALSE, quote = FALSE)

#report <- read.csv(paste0(cluster_dir, "raw_report.csv")

# Remove controls

control_report <- report %>% 
  dplyr::left_join(controls, by = c("run", "sample")) %>%
  dplyr::filter(is.na(control)) %>% # just keep the controls
  dplyr::select(-control, dnarna) %>% 
  dplyr::group_by(sample, run, db, species, species_taxid, type, tool) %>% 
  dplyr::summarise(control_reads = sum(reads)) %>% 
  dplyr::left_join(read_counts, by = c("sample", "run")) %>% 
  dplyr::mutate(control_rpm = control_reads*1000000/raw_reads) %>%
  dplyr::rename(control = sample, control_raw_reads = raw_reads)

report_controls_rm <- report %>% 
  dplyr::left_join(controls, by = c("run", "sample")) %>% 
  dplyr::filter(!is.na(control)) %>% 
  dplyr::rename(sample_reads = reads) %>% 
  dplyr::left_join(read_counts, by = c("sample", "run")) %>% 
  dplyr::mutate(sample_rpm = sample_reads*1000000/raw_reads) %>% 
  dplyr::rename(sample_raw_reads = raw_reads) %>% 
  dplyr::left_join(control_report, by = c("run", "db", "species", "species_taxid", "type", "tool", "control")) %>% 
  dplyr::mutate(control_rpm = ifelse(is.na(control_reads) | control_reads == 0, 1000000/sample_raw_reads, control_rpm)) %>%  # set control rpm to equivalent of 1 read in the sample 
  dplyr::mutate(rpm_ratio = sample_rpm / control_rpm)

# Thresholds

# Remove taxids in human lineage
human_taxids <- c(131567, 2759, 33154, 33208, 6072, 33213, 33511, 7711, 89593, 7742, 7776, 117570, 117571,
                   8287, 1338369, 32523, 32524, 40674, 32525, 9347, 1437010, 314146, 9443, 376913, 9526, 314295, 9604, 207598, 9605, 9606)

human_taxons <- c("cellular organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Deuterostomia", "Chordata", 
"Craniata", "Vertebrata", "Gnathostomata", "Teleostomi", "Euteleostomi", "Sarcopterygii", "Dipnotetrapodomorpha", "Tetrapoda",
"Amniota", "Mammalia", "Theria", "Eutheria", "Boreoeutheria", "Euarchontoglires", "Primates", "Haplorrhini", "Simiiformes", "Catarrhini",
"Hominoidea", "Hominidae", "Homininae", "Homo", "Homo sapiens")

totals_nonhuman_classified <- report_controls_rm %>% 
  dplyr::filter(!(species %in% c("unclassified", "root", "unidentified", "uncultured microorganism", "uncultured organism"))
                & !(species %in% human_taxons)) %>% 
  dplyr::group_by(sample, run, db, tool) %>% 
  dplyr::summarise(total_sample_reads_nonhuman_classified = sum(sample_reads),
                   total_control_reads_nonhuman_classified = sum(control_reads, na.rm = TRUE))

totals_species <- report_controls_rm %>% 
  dplyr::filter(!is.na(type) & type != "Archaea" & !is.na(species_taxid)) %>% 
  dplyr::group_by(sample, run, db, tool) %>% 
  dplyr::summarise(total_sample_reads_species = sum(sample_reads))

totals_type <- report_controls_rm %>% 
  dplyr::group_by(sample, run, db, tool, type) %>% 
  dplyr::summarise(total_sample_reads_type = sum(sample_reads))

report_clean <- report_controls_rm %>% 
  dplyr::left_join(totals_nonhuman_classified, by = c("sample", "run", "db", "tool")) %>% 
  dplyr::left_join(totals_species, by = c("sample", "run", "db", "tool")) %>% 
  dplyr::left_join(totals_type, by = c("sample", "run", "db", "tool", "type")) %>% 
  dplyr::mutate(proportion_nonhuman_classified = sample_reads / total_sample_reads_nonhuman_classified,
                proportion_species = sample_reads / total_sample_reads_species,
                proportion_type = sample_reads / total_sample_reads_type) %>% 
  dplyr::mutate(control_reads = coalesce(control_reads, 0)) %>% 
  dplyr::filter(sample_reads != 0)

write.csv(report_clean, paste0(cluster_dir, "full_report.csv"), row.names = FALSE)

# Just species in mock community
msa_2008_species <- community_taxid_mapping %>% 
  dplyr::filter(community == "msa_2008")

report_mc <- report_clean %>% 
  dplyr::filter(species_taxid %in% msa_2008_species$taxid | species_taxid %in% msa_2008_species$taxid_czid)

write.csv(report_mc, paste0(cluster_dir, "mock_community_report.csv"), row.names = FALSE, quote = FALSE)
  

