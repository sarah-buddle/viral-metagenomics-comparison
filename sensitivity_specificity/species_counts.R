library(tidyverse)
library(taxonomizr)

options(scipen=999)

# Thresholds to test
identifiers_list <- c(5)

# Inputs
report_clean <- read.csv(report_input) %>% 
  dplyr::mutate(db = replace_na(db, "NA")) %>% 
  dplyr::filter(!is.na(species_taxid))

thresholds <-  read.csv(thresholds_input)


# Mappings  
community_taxid_mapping <- read.csv(paste0(local_dir, "community_taxid_mapping.csv"))

run_community_mapping <- read.csv(paste0(local_dir, "run_community_mapping.csv"))

mappings <- full_join(community_taxid_mapping, run_community_mapping, by = "community") %>% 
  dplyr::rename(species = name)

controls <- read.csv(paste0(local_dir, "corresponding_controls.csv"))

# Function to determine presence/absence

isPositive <- function(type, rpm_ratio, sample_reads, control_reads, proportion, thres_ref_no) {
  
  thres_list <- thresholds %>% 
    dplyr::filter(thres_ref == thres_ref_no)
  
  rpm_ratio_thres_bac <- as.numeric(thres_list$rpm_ratio_thres_bac[[1]])
  rpm_ratio_thres_vir <- as.numeric(thres_list$rpm_ratio_thres_vir[[1]])
  rpm_ratio_thres_vir_higher <- as.numeric(thres_list$rpm_ratio_thres_vir_higher[[1]])
  rpm_ratio_thres_fun <- as.numeric(thres_list$rpm_ratio_thres_fun[[1]])
  rpm_ratio_thres_euk <- as.numeric(thres_list$rpm_ratio_thres_euk[[1]])
  proportion_thres_bac <- as.numeric(thres_list$proportion_thres_bac[[1]])
  proportion_thres_vir <- as.numeric(thres_list$proportion_thres_vir[[1]])
  proportion_thres_fun <- as.numeric(thres_list$proportion_thres_fun[[1]])
  proportion_thres_euk <- as.numeric(thres_list$proportion_thres_euk[[1]])
  reads_thres_vir <- as.numeric(thres_list$reads_thres_vir[[1]])
  reads_thres_fun <-  as.numeric(thres_list$reads_thres_vir[[1]])
  control_reads_thres_vir <- as.numeric(thres_list$control_reads_thres_vir[[1]])
  
  output <- NA
  
  if (!is.na(type) & !is.na(rpm_ratio) & !is.na(sample_reads) & !is.na(proportion)) {
  
    if (type == "Bacteria"){
      
      if (rpm_ratio >= rpm_ratio_thres_bac & proportion >= proportion_thres_bac) {
        
        output <- "positive"
        
      } else {
        
        output <- "negative"
        
      }
      
    } else if (type == "Virus") {
      
      if (rpm_ratio >= rpm_ratio_thres_vir & proportion >= proportion_thres_vir &
          sample_reads >= reads_thres_vir & 
          (control_reads <= control_reads_thres_vir | rpm_ratio >= rpm_ratio_thres_vir_higher )) {
        
        output <- "positive"
        
      } else {
        
        output <- "negative"
        
      }
      
    } else if (type == "Fungi") {
      
      if (rpm_ratio >= rpm_ratio_thres_fun & sample_reads >= reads_thres_fun &
          proportion >= proportion_thres_fun) {
        
        output <- "positive"
        
      } else {
        
        output <- "negative"
        
      }
      
    } else if (type == "Other eukaryote") {
      
      if (rpm_ratio >= rpm_ratio_thres_euk & proportion >= proportion_thres_euk) {
        
        output <- "positive"
        
      } else {
        
        output <- "negative"
        
      }
      
    }
    
  } else {
    
    output <- NA
    
  }
  
  return(output)
  
}

isPositive <- Vectorize(isPositive, vectorize.args = c("type", "rpm_ratio", "sample_reads", "control_reads",
                                                       "proportion", "thres_ref_no"))

# Dataframe of all positve species
positive_species <- report_clean %>% 
  dplyr::filter(!(tool %in% c("minimap", "minimap_nodedup", "bowtie", "bowtie_nodedup"))) %>% 
  dplyr::mutate(thres_ref_no = map(row_number(), ~ rep(identifiers_list, each = 1))) %>%
  tidyr::unnest(thres_ref_no) %>% 
  dplyr::mutate(result = isPositive(type, rpm_ratio, sample_reads, control_reads, proportion_nonhuman_classified, thres_ref_no)) %>% # use the type of proportion needed here
  dplyr::filter(result == "positive")

positive_species_stats <- read.csv(positive_species_stats_input)

positive_species_nothres <- report_clean %>% 
  dplyr::filter(!(tool %in% c("minimap", "bowtie"))) %>% 
  dplyr::mutate(thres_ref_no = 0,
                result = "positive")

positive_species <- rbind(positive_species, positive_species_stats, positive_species_nothres)

# Overwrite
write.csv(positive_species, positive_species_output, row.names = FALSE, quote = FALSE)

# Combine results from DNA and RNA where separate
positive_species_dnarna <- positive_species %>% 
  arrange(proportion_nonhuman_classified) %>% 
  distinct(run, dnarna_pair, species, species_taxid, type, db, tool, thres_ref_no, repeat., dilution)

# False by type e.g. bacteria, viruses
msa_2008 <- community_taxid_mapping %>% 
  filter(community == "msa_2008")

expected_taxids <- unique(c(msa_2008$taxid, msa_2008$taxid_czid))

false_positives_bytype <- positive_species_dnarna %>%
  dplyr::mutate(count = 1) %>%
  dplyr::filter(!(species_taxid %in% expected_taxids)) %>% 
  dplyr::group_by(dnarna_pair, run, db, repeat., dilution, type, tool, thres_ref_no) %>% 
  dplyr::summarise(false_positives = sum(count))

write.csv(false_positives_bytype, false_positives_bytype_output, row.names = FALSE, quote = FALSE)

# Count species
species_counts <- positive_species_dnarna %>%
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(run, dnarna_pair, db, repeat., dilution, tool, thres_ref_no) %>% 
  dplyr::summarise(total_species = sum(count)) %>% 
  dplyr::mutate(type = "All")
  
# True positive counts
true_positives <- data.frame()

for (run_name in unique(positive_species_dnarna$run)) {
  
  # Get expected taxids for that run
  run_mappings <- mappings %>% 
    dplyr::filter(run == run_name)
  
  run_mappings_mock <- run_mappings %>% 
    filter(positive_control == FALSE)
  
  run_mappings_control <- run_mappings %>% 
    filter(positive_control == TRUE)
  
  # Filter to just relevant records
  
  run_report <- positive_species_dnarna %>% 
    dplyr::filter(run == run_name)
  
  for (tool_name in unique(run_report$tool)) {
    
    tool_report <- run_report %>% 
      dplyr::filter(tool == tool_name)
    
    expected_taxids <- NA
    
    if (tool_name %in% c("czid", "czid_nt_count", "czid_nr_count")) {
      
      expected_taxids_mock <- run_mappings_mock$taxid_czid
      
      expected_taxids_control <- run_mappings_control$taxid_czid
      
    } else {
      
      expected_taxids_mock <- run_mappings_mock$taxid
      
      expected_taxids_control <- run_mappings_control$taxid
      
    }
  
    for (sample_name in unique(tool_report$dnarna_pair)) {
      
      sample_report <- tool_report %>% 
        dplyr::filter(dnarna_pair == sample_name)
      
      for (db_name in unique(sample_report$db)){
        
        db_report <- sample_report %>% 
          dplyr::filter(db == db_name)
        
        for (thres_ref_no_name in unique(db_report$thres_ref_no)) {
          
          thres_report <- db_report %>% 
            dplyr::filter(thres_ref_no == thres_ref_no_name)
          
          # Count true positives
          all_taxids <- unique(thres_report$species_taxid)
      
          true_pos <- 0
          true_pos_taxids <- c()
          for (taxid in expected_taxids_mock) {
            if (taxid %in% all_taxids) {
              true_pos <- true_pos + 1
              true_pos_taxids <- c(true_pos_taxids, taxid)
            } 
          }
          
          # Count controls
          pos_controls <- 0
          pos_control_taxids <- c()
          for (taxid in expected_taxids_control) {
            if (taxid %in% all_taxids) {
              pos_controls <- pos_controls + 1
              pos_control_taxids <- c(pos_control_taxids, taxid)
            } 
            
          }
          
          # Merge with species counts data frame
          new_row <- data.frame(sample_name, run_name, db_name, tool_name, thres_ref_no_name, true_pos, pos_controls)
          new_row$true_pos_taxids <- paste(true_pos_taxids, collapse = " ")
          new_row$pos_control_taxids <- paste(pos_control_taxids, collapse = " ")
          new_row$pos_taxids <- paste(all_taxids, collapse = " ")
          names(new_row) <- c("dnarna_pair", "run", "db", "tool", "thres_ref_no", "true_positives", "positive_controls",
                              "true_positive_taxids", "positive_control_taxids", "positive_taxids")
          
          true_positives <- rbind(true_positives, new_row)
          
        }
      }
    }
  }
}

species_counts_final <- species_counts %>% 
  dplyr::full_join(true_positives, by = c("dnarna_pair", "run", "db", "tool", "thres_ref_no")) %>% 
  mutate(t_pos = trimws(paste(positive_control_taxids, true_positive_taxids, sep = " "))) %>%
  mutate(t_pos = trimws(gsub(" ", "|", t_pos)))

# Find false positive taxids
species_counts_final$false_positive_taxids <- trimws(sapply(1:nrow(species_counts_final), 
                                                             function(x) gsub(species_counts_final$t_pos[x], "", 
                                                             species_counts_final$positive_taxids[x])))

theoretical <- function(community) {
  
  if (community == "msa_2008") {
    return(6)
  } else if (community == "zymo_d6306") {
    return(10)
  } else {
    return(NA)
  }
  
}

theoretical <- Vectorize(theoretical)

find_names <- function(taxids) {
  
  names <- c()
  
  if (taxids != "") {
  
    positive_taxids_list <- as.list(strsplit(taxids, " "))[[1]]
    
    for (taxid in positive_taxids_list) {
      
      if (taxid != "") {
        
        name_list <- getCommon(taxid, types = "scientific name")
        
        name <- name_list[[1]][[1]]
        
        names <- c(names, name)
        
      }
      
    }
    
    return(paste(names, collapse = " "))
    
  } else {
    
    return("")
    
  }
  
}

find_names <- Vectorize(find_names)

species_counts_final2 <- species_counts_final %>%
  select(-t_pos) %>% 
  left_join(run_community_mapping, by = "run") %>% 
  mutate(theoretical = theoretical(community),
         false_positives = total_species - true_positives - positive_controls,
         false_negatives = theoretical - true_positives,
         precision = true_positives / (true_positives + false_positives),
         recall = true_positives / (true_positives + false_negatives),
         fscore = (2*precision*recall)/(precision+recall)) %>% 
  mutate(false_negatives = ifelse(false_negatives < 0, 0, false_negatives))

write.csv(species_counts_final2, species_counts_output, row.names = FALSE, quote = FALSE)



