#### Stats for thresholds ####

library(tidyverse)
library(tidymodels)
library(caret)

runs <- c("illumina_250923", "nanopore_270923", "twist_101123")

tool_levels_main <- c("kraken", "bracken", "dragen", "epi2me_bracken", "megan_nucleotide",
                      "metamix_nucleotide", "metamix_step2_nucleotide", "czid_nt_count", "onecodex")

dilutions <- c(10, 100, 1000, 10000)

## File paths ##
positive_species_stats_output <- paste0(cluster_dir, "positive_species_stats.csv")

community_taxid_mapping <- read.csv(paste0(local_dir, "community_taxid_mapping.csv")

report_input <- paste0(cluster_dir, "full_report.csv")

## Import and clean report from full_report ##
report_raw <- read.csv(report_input)

expected_taxids <- community_taxid_mapping %>% 
  filter(community == "msa_2008")

true_positives <- unique(c(expected_taxids$taxid, expected_taxids$taxid_czid))

phage_taxids <- c(329852, 2681611, 10710)
                  

is_community <- function(species_taxid, expected_taxids) {
  
  if (species_taxid %in% expected_taxids) {
    return("positive")
  } else {
    return("negative")
  }
 
}

is_community <- Vectorize(is_community, vectorize.args = "species_taxid")

report <- report_raw %>% 
  dplyr::filter(type == "Virus") %>% 
  dplyr::filter(tool %in% tool_levels_main & grepl("_sub", sample) & run %in% runs & dilution %in% dilutions &
  !(species_taxid %in% phage_taxids) & !is.na(species_taxid)) %>% 
  dplyr::mutate(expected = as.factor(is_community(species_taxid, true_positives))) %>% 
  dplyr::filter(!(run == "twist_101123" & tool == "onecodex"))

#### Individual ROC curves - SI Figure 2 ####

# ROC curve for rpmr only

roc_rpmr <- report %>% 
  tidymodels::roc_curve(., truth = expected, rpm_ratio, event_level = "second") %>% 
  dplyr::mutate(difference = sensitivity - specificity) %>% 
  dplyr::mutate(name = "rpmr_only") %>% 
  dplyr::rename(threshold = .threshold)

my_roc_auc_rpmr <- roc_auc(report, truth = expected, rpm_ratio, event_level = "second")

auc_rpmr <- my_roc_auc_rpmr$.estimate[1] %>% 
  round(., digits = 3)

# Find threshold where sensitivity = specificity
rpmr_threshold_df_ss <- roc_rpmr %>%
  dplyr::filter(difference > 0) %>% 
  dplyr::arrange(difference) %>% 
  head(., 1)

rpmr_threshold_ss <- rpmr_threshold_df_ss$threshold[1]

# Find the threshold where sensitivity = 95%
rpmr_threshold_95_df <- roc_rpmr %>% 
  filter(sensitivity >= 0.95) %>% 
  arrange(sensitivity) %>% 
  head(., 1)

rpmr_threshold_95 <- rpmr_threshold_95_df$threshold[1]

# Use 5 as threshold
rpmr_threshold_5_df <- roc_rpmr %>% 
  filter(threshold >= 5) %>% 
  arrange(threshold) %>% 
  head(., 1)

rpmr_threshold_5 <- rpmr_threshold_5_df$threshold[1]

# Function to calculate roc when some previous results have been excluded
my_roc <- function(report, truth, variable1, variable2, variable1_thres) {
  
  report_pass <- report %>% 
    dplyr::filter(is.na(control_reads) | control_reads == 0 | rpm_ratio >= variable1_thres)
  
  report_reject <- report %>% 
    dplyr::filter(!(is.na(control_reads) | control_reads == 0 | rpm_ratio >= variable1_thres))
  
  report_reject_n <- as.numeric(nrow(report_reject))
  
  previous_fn <- report_reject %>% 
    dplyr::filter(expected == "positive") %>% 
    nrow(.) %>% 
    as.numeric(.)
  
  previous_tn <- report_reject_n - previous_fn
  
  df <- report_pass %>% 
    dplyr::select(all_of(c(truth, variable2))) %>% 
    dplyr::arrange(!!as.symbol(variable2))
  
  names <- c("threshold", "sensitivity", "specificity")
  
  output <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(output) <- names
  
  for (threshold in c(0, df[[variable2]])) {
    
    df_pred <- df %>% 
      dplyr::mutate(pred = ifelse(!!as.symbol(variable2) >= threshold, "positive", "negative")) %>% 
      dplyr::mutate(truth = factor(!!as.symbol(truth), c("positive", "negative")),
             pred = factor(pred, c("positive", "negative")))
    
    cm <- conf_mat(df_pred, truth = truth, estimate = pred)
    
    tp <- cm$table[1,1]
    fn <- cm$table[2,1]
    fp <- cm$table[1,2]
    tn <- cm$table[2,2]
    
    sensitivity <- tp / (tp + fn + previous_fn)
    specificity <-  (tn + previous_tn) / (fp + tn + previous_tn)
    
    new_row <- data.frame(t(data.frame(c(threshold, sensitivity, specificity))))
    
    colnames(new_row) <- names
    rownames(new_row) <-  NULL
    
    output <- bind_rows(output, new_row)
    
  }
  
  output_final <- output %>% 
    dplyr::mutate(difference = sensitivity - specificity)
  
  return(output_final)
  
}

#### Sensitivity = Speciicity ####

roc_pmr_ss <- my_roc(report, "expected", "rpm_ratio", "proportion_nonhuman_classified", rpmr_threshold_ss) %>%  
  dplyr::mutate(name = "rpmr+pmr_s=s")

pmr_threshold_df_ss <- roc_pmr_ss %>%
  dplyr::filter(difference > 0) %>% 
  dplyr::arrange(difference) %>% 
  head(., 1)

pmr_threshold_ss <- pmr_threshold_df_ss$.threshold[1]


#### Sensitivity = 95% ####
roc_pmr_95 <- my_roc(report, "expected", "rpm_ratio", "proportion_nonhuman_classified", rpmr_threshold_95) %>%  
  dplyr::mutate(name = "rpmr+pmr_95")

pmr_threshold_df_95 <- roc_pmr_95 %>%
  dplyr::filter(sensitivity >= 0.95) %>% 
  dplyr::arrange(sensitivity) %>% 
  head(., 1)

pmr_threshold_95 <- pmr_threshold_df_95$.threshold[1]

# Apply our thresholds
roc_pmr_5 <- my_roc(report, "expected", "rpm_ratio", "proportion_nonhuman_classified", rpmr_threshold_5) %>%  
  dplyr::mutate(name = "rpmr=5+pmr=0.0001")

pmr_threshold_df_5 <- roc_pmr_5 %>%
  dplyr::filter(threshold >= 0.0001) %>% 
  dplyr::arrange(threshold) %>% 
  head(., 1)

pmr_threshold_5 <- pmr_threshold_df_5$.threshold[1]

## Fit logistic regression model ##
model <- tidymodels::logistic_reg() %>% 
  tidymodels::set_engine("glm") %>% 
  tidymodels::fit(expected ~ proportion_nonhuman_classified + rpm_ratio, data = report)

preds <- model %>% 
  tidymodels::augment(new_data = report)

roc_lm <- tidymodels::roc_curve(preds, truth = expected, .pred_positive, event_level = "second") %>% 
  dplyr::mutate(difference = sensitivity - specificity,
         name = "rpmr+pmr_lm") %>% 
  dplyr::rename(threshold = .threshold)

# sensitivity = specificity
lm_threshold_df_ss <- roc_lm %>%
  dplyr::filter(difference >= 0) %>% 
  dplyr::arrange(difference) %>% 
  head(., 1)

lm_threshold_ss <- lm_threshold_df_ss$threshold[1]

# sensitivity > 0.95
lm_threshold_df_95 <- roc_lm %>%
  dplyr::filter(sensitivity >= 0.95) %>% 
  dplyr::arrange(sensitivity) %>% 
  head(., 1)

lm_threshold_95 <- lm_threshold_df_95$threshold[1]

# AUC
my_roc_auc <- roc_auc(preds, truth = expected, .pred_positive, event_level = "second")

auc <- my_roc_auc$.estimate[1] %>% 
  round(., digits = 3)

## Import decontam ROC curve ##
roc_decontam <- read.csv(paste0(cluster_dir, "decontam_roc.csv") %>% 
  dplyr::select(threshold, sensitivity, specificity) %>% 
  dplyr::mutate(difference = sensitivity - specificity,
         name = "decontam") %>% 
  dplyr::distinct()

name_levels <- c("rpmr_only", "rpmr=5+pmr=0.0001", "rpmr+pmr_lm", "decontam")

name_labels <- c("Reads per million ratio only", "Read per million ratio and\nproportion of microbial reads",
                  "Linear model", "Decontam")

## Plot ROC curves ##
rbind(roc_rpmr, roc_pmr_5, roc_lm, roc_decontam) %>% 
  dplyr::mutate(name = factor(name, name_levels, name_labels)) %>% 
  ggplot2::ggplot(aes(x = 1 - specificity, y = sensitivity, color = name)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = "right",
        legend.text = element_text(size = 8))

ggsave(paste0(local_dir, "/virus_methods/roc_curves/roc_curve_combined_decontam.png"),
       width = 5, height = 3)

rbind(roc_rpmr, roc_pmr_5, roc_lm) %>% 
  mutate(name = factor(name, name_levels, name_labels)) %>% 
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = name)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.position = "right",
        legend.text = element_text(size = 8))

ggsave(paste0(local_dir, "virus_methods/roc_curves/roc_curve_combined.png",
       width = 5, height = 3)


#### RPMR only  #####
rpmr_thres5 <- report %>% 
  dplyr::mutate(estimate = ifelse(rpm_ratio >= 5 | is.na(control_reads) | control_reads == 0, "positive", "negative")) %>% 
  dplyr::mutate(estimate = as.factor(estimate))

cm <- caret::confusionMatrix(data = rpmr_thres5$estimate, reference = rpmr_thres5$expected, positive = "positive")

cm$byClass["Sensitivity"]

cm$byClass["Specificity"]



  

