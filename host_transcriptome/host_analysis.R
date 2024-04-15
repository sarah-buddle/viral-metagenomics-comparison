library(tidyverse)
library(tximport)
library(rtracklayer)


dir.create(outdir)
dir.create(paste0(outdir, "/gene_lists"))

# Import human annotation - basic version doesn't include some transcript IDs
gtf <- rtracklayer::import(human_gtf)

gtf_df <- as.data.frame(gtf) %>% 
  dplyr::filter(type == "gene")
 
# Link transcript ID to gene ID
tx2gene <- as.data.frame(gtf) %>% 
  dplyr::select(transcript_id, gene_id) %>% 
  dplyr::distinct() %>% 
  dplyr::filter(!is.na(transcript_id))

#### Import kallisto outputs ####

import_kallisto <- function(samples, platform, run, kallisto_filepath, tx2gene) {
  
  cpm <- data.frame()
  
  for (sample in samples) {
    
    input <- paste0(kallisto_filepath, "/", run, "/", sample, "/abundance_newname.tsv")
    
    names(input) <- "counts"
    
    kallisto <- tximport(input, type = "kallisto", tx2gene = tx2gene, countsFromAbundance ="no")
    
    kallisto_df <- as.data.frame(kallisto$counts)
    colnames(kallisto_df) <- colnames(kallisto$counts)
    
    total <- sum(kallisto$counts)
    
    cpm_sample <- kallisto_df %>% 
      tibble::rownames_to_column(var = "feature") %>% 
      dplyr::mutate(cpm = counts*10^6 / total,
             sample = sample,
             run = run,
             platform = platform)
    
    cpm <- rbind(cpm, cpm_sample)
    
  }
  
  cpm_mean <- cpm %>% 
    dplyr::group_by(feature, platform) %>% 
    dplyr::summarise(mean_cpm = mean(cpm),
              mean_counts = mean(counts))
  
  
}

illumina_samples <- c("RNA_10000_a_sub", "RNA_10000_b_sub", "RNA_1000_a_sub", "RNA_1000_b_sub",
                      "RNA_100_a_sub", "RNA_100_b_sub", "RNA_10_a_sub", "RNA_10_b_sub")

kallisto_illumina <- import_kallisto(illumina_samples, "illumina", "illumina_250923", kallisto_dir, tx2gene)

nanopore_samples <- c("barcode11_sub", "barcode12_sub", "barcode13_sub", "barcode14_sub", "barcode15_sub", "barcode16_sub",
                      "barcode17_sub", "barcode18_sub")

kallisto_nanopore <- import_kallisto(nanopore_samples, "nanopore", "nanopore_270923", kallisto_dir, tx2gene)

twist_samples <- c("10000_a_sub", "10000_b_sub", "1000_a_sub", "1000_b_sub",
                   "100_a_sub", "100_b_sub", "100_c_sub", "100_d_sub",
                   "10000_a_sub", "10_b_sub", "10_c_sub", "10_d_sub")

kallisto_twist <- import_kallisto(twist_samples, "twist", "twist_101123", kallisto_dir, tx2gene)

coding_levels <- c("Protein-coding", "Non-coding")

# Add gene type and biotype data from gtf
expression <- rbind(kallisto_illumina, kallisto_nanopore, kallisto_twist) %>%
  dplyr::left_join(gtf_df, join_by(feature == gene_id)) %>% 
  dplyr::select(feature, platform, mean_cpm, mean_counts, type, gene_type) %>% 
  dplyr::filter(type == "gene") %>% 
  dplyr::mutate(coding = ifelse(gene_type == "protein_coding", "Protein-coding", "Non-coding")) %>% 
  dplyr::mutate(coding = factor(coding, coding_levels)) %>% 
  dplyr::mutate(log2_cpm = log2(mean_cpm))

write.csv(expression, paste0(outdir, "/expression.csv"), 
          quote = FALSE, row.names = FALSE)

platform_levels <- c("illumina", "nanopore", "twist")

platform_labels <- c("Illumina", "ONT", "Twist_VRP")

expression <- read.csv(paste0(outdir, "/expression.csv")) %>% 
  dplyr::mutate(coding = factor(coding, coding_levels),
                platform = factor(platform, platform_levels, platform_labels))

#### CPM scatterplots ####

cpm_scatterplot <- function(expression, platform_x, platform_y, outdir) {
  
  expression_xy <- expression %>% 
    dplyr::select(-c(mean_cpm, mean_counts)) %>% 
    tidyr::pivot_wider(names_from = platform, values_from = log2_cpm, values_fill = -Inf) %>% 
    dplyr::filter(!!as.symbol(platform_x) != -Inf | !!as.symbol(platform_y) != -Inf)
  
  p_cor <- cor.test(expression_xy[[platform_x]], expression_xy[[platform_y]], 
                    method = "spearman")$estimate %>% 
    round(., digits = 3)
  
  expression_xy_pc <- expression_xy %>% 
    dplyr::filter(coding == "Protein-coding")
  
  p_cor_pc <- cor.test(expression_xy_pc[[platform_x]], expression_xy_pc[[platform_y]], 
                    method = "spearman")$estimate %>% 
    round(., digits = 3)
  
  plot <- ggplot2::ggplot(expression_xy, aes(x = !!as.symbol(platform_x), y = !!as.symbol(platform_y), color = coding)) +
    geom_point(alpha = 0.1, size = 0.05) +
    theme(panel.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "top",
          plot.caption = element_text(hjust = 0.5)) +
    scale_x_continuous(expand = expansion(mult = c(0, .1)), limits = c(-10, 15)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(-10, 15)) +
    xlab(paste0(str_to_title(gsub("_", " ", platform_x)), " log2(CPM)")) +
    ylab(paste0(str_to_title(gsub("_", " ", platform_y)), " log2(CPM)")) +
    geom_hline(yintercept = 0, linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(color = "Gene biotype") +
    guides(colour = guide_legend(override.aes = list(size =3 , alpha = 1))) +
    labs(caption = paste0("Spearman rho: ", p_cor, " (all)      ", p_cor_pc, " (protein coding only)"))
  
  output_filepath <- paste0(outdir, "/", platform_x, "_", platform_y, "_cpm.png")
  
  ggsave(output_filepath, plot, height = 4.2, width = 4)
  
}

cpm_scatterplot(expression, "Illumina", "ONT", outdir)

cpm_scatterplot(expression, "Illumina", "Twist_VRP", outdir)

cpm_scatterplot(expression, "ONT", "Twist_VRP", outdir)


#### Total counts ####

total_counts <- expression %>% 
  dplyr::group_by(platform) %>% 
  dplyr::summarise(total_mean_counts = sum(mean_counts))

write.csv(total_counts, paste0(outdir, "/total_counts.csv"), quote = FALSE, row.names = FALSE)

#### Bar plot ####

read_thresholds <- function(mean_counts, platform, threshold_illumina, threshold_nanopore) {
  
  if (platform %in% c("Illumina", "Twist_VRP") & mean_counts >= threshold_illumina) {
    
    return(mean_counts)
    
  } else if (platform == "ONT" & mean_counts >= threshold_nanopore) {
    
    return(mean_counts)
    
  } else {
    
    return(0)
    
  }
  
}

read_thresholds <- Vectorize(read_thresholds,vectorize.args = c("mean_counts", "platform"))

platform_levels <- c("Illumina", "ONT", "Twist_VRP", "Both")

platform_labels <- c("Illumina only", "ONT only", "Twist Panel only", "Both")

find_platform <- function(x, y, x_name, y_name) {
  
  if (!is.na(x) & !is.na(y) & x != 0 & y != 0) {
    
    return("Both")
    
  } else if  (!is.na(x) & x != 0 ) {
    
    return(x_name)
    
  } else if  (!is.na(y) & y != 0 ) {
    
    return(y_name)
    
  }
  
}

find_platform <- Vectorize(find_platform, vectorize.args = c("x", "y"))

gene_list <- function(coding_counts, platform1, platform2, found_by_name, coding_name, threshold_illumina, threshold_nanopore) {
  
  genes_df <- coding_counts %>% 
    dplyr::filter(found_by == found_by_name & coding == coding_name)
  
  genes <- unique(genes_df$feature)
  
  filepath <- paste0(outdir, "/gene_lists/", platform1, "_", platform2, "_", found_by_name, "_", coding_name, "_ti", threshold_illumina, "_tn", threshold_nanopore, ".txt")
  
  file.remove(filepath)
  
  lapply(genes, write, filepath, append = TRUE, ncolumns = 1)
  
}

coding_barplot <- function(expression, platform1, platform2, threshold_illumina, threshold_nanopore, outdir, platform_colours, ylim) {
  
  coding_counts <- expression %>% 
    dplyr::select(-c(mean_cpm, log2_cpm)) %>% 
    dplyr::mutate(mean_counts = read_thresholds(mean_counts, platform, threshold_illumina, threshold_nanopore)) %>% 
    dplyr::filter(platform %in% c(platform1, platform2)) %>% 
    tidyr::pivot_wider(names_from = "platform", values_from = "mean_counts") %>% 
    dplyr::filter(!!as.symbol(platform1) != 0 | !!as.symbol(platform2) != 0) %>% 
    dplyr::mutate(found_by = find_platform(!!as.symbol(platform1), !!as.symbol(platform2), platform1, platform2),
           count = 1) 
  
  # Make bar plot
  bar <- coding_counts %>% 
    dplyr::group_by(coding, found_by) %>% 
    dplyr::summarise(n_genes = sum(count)) %>% 
    dplyr::mutate(found_by = factor(found_by, platform_levels, platform_labels))
  
  
  bar_output <- paste0(outdir, "/", platform1, "_", platform2, "_ti", threshold_illumina, "_tn", threshold_nanopore, "_coding.png")
  
  bar %>%   
    ggplot(aes(x = coding, y = n_genes, fill = found_by)) +
    geom_bar(stat = "identity") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0, ylim)) +
    ylab("Number of genes")+
    scale_fill_manual(values = platform_colours)
  
  ggsave(bar_output, height = 3, width = 1.5)
  
  # Make gene lists for GO analysis
  
 gene_list(coding_counts, platform1, platform2, "Both", "Protein-coding", threshold_illumina, threshold_nanopore)
 gene_list(coding_counts, platform1, platform2, "Both", "Non-coding", threshold_illumina, threshold_nanopore)
 gene_list(coding_counts, platform1, platform2, platform1, "Protein-coding", threshold_illumina, threshold_nanopore)
 gene_list(coding_counts, platform1, platform2, platform1, "Non-coding", threshold_illumina, threshold_nanopore)
 gene_list(coding_counts, platform1, platform2, platform2, "Protein-coding", threshold_illumina, threshold_nanopore)
 gene_list(coding_counts, platform1, platform2, platform2, "Non-coding", threshold_illumina, threshold_nanopore)
  
}

platform_colours2 <- setNames(c("#D95F02", "#7570B3", "#E7298A", "#1B9E77"),
                              c("Illumina only", "ONT only", "Twist Panel only", "Both"))

coding_barplot(expression, "Illumina", "ONT", threshold_illumina=0, threshold_nanopore=0, outdir, platform_colours2, 40000)
coding_barplot(expression, "Illumina", "ONT", threshold_illumina=10, threshold_nanopore=3, outdir, platform_colours2, 20000)
coding_barplot(expression, "Illumina", "Twist_VRP", threshold_illumina=0, threshold_nanopore=0, outdir, platform_colours2, 40000)
coding_barplot(expression, "Illumina", "Twist_VRP", threshold_illumina=10, threshold_nanopore=3, outdir, platform_colours2, 20000)
coding_barplot(expression, "ONT", "Twist_VRP", threshold_illumina=0, threshold_nanopore=0, outdir, platform_colours2, 40000)
coding_barplot(expression, "ONT", "Twist_VRP", threshold_illumina=10, threshold_nanopore=3, outdir, platform_colours2, 20000)


#### Plot showing CPM and how many platforms each point is detected by ####

which_platforms <- function(illumina, nanopore, twist) {
  
  if (illumina != -Inf & nanopore != -Inf & twist != -Inf) {
    
    return("All")
    
  } else if (illumina != -Inf & nanopore != -Inf) {
    
    return("Illumina + ONT")
    
  } else if (illumina != -Inf & twist != -Inf) {
    
    return("Illumina + Twist VRP")
    
  } else if (nanopore != -Inf & twist != -Inf) {
    
    return("ONT + Twist VRP")
    
  } else if (illumina != -Inf) {
    
    return("Illumina only")
    
  } else if (nanopore != -Inf) {
    
    return("ONT only")
    
  } else if (twist != -Inf) {
    
    return("Twist VRP only")
    
  } else {
    
    return("None")
    
  }
  
}

which_platforms <- Vectorize(which_platforms)

platform_comb_labels <- c("All", "Illumina + ONT", "Illumina + Twist VRP", "ONT + Twist VRP",
                          "Illumina only", "ONT only", "Twist VRP only", "None")

platform_levels <- c("Illumina", "ONT", "Twist_VRP")

platform_labels <- c("Illumina log2(CPM)", "ONT log2(CPM)", "Twist VRP log2(CPM)")


expression_platform <- expression %>% 
  dplyr::select(-c(mean_cpm, mean_counts)) %>% 
  tidyr::pivot_wider(names_from = "platform", values_from = "log2_cpm") %>% 
  dplyr::mutate(which_platforms = which_platforms(Illumina, ONT, Twist_VRP)) %>% 
  dplyr::filter(which_platforms != "None") %>% 
  tidyr::pivot_longer(cols = c(Illumina, ONT, Twist_VRP), names_to = "platform", values_to = "log2_cpm") %>% 
  dplyr::mutate(which_platforms = factor(which_platforms, levels = platform_comb_labels, labels = platform_comb_labels),
         platform = factor(platform, platform_levels, platform_labels))


expression_platform %>% 
  ggplot2::ggplot(aes(x = which_platforms, y = log2_cpm, color = which_platforms)) +
  geom_boxplot(outlier.color = "white") +
  facet_grid(cols = vars(platform), space = "free_x", scales = "free_x") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.key = element_blank(),
        legend.position = "right",
        plot.caption = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),
        legend.title = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(-12, 12)) +
  ylab("log2(CPM)") +
  scale_color_brewer(palette = "Dark2") +
  labs(caption = "outliers not shown")

ggsave(paste0(outdir, "/platform_boxplot.png"), height = 3.7, width = 9, dpi = 1000)

expression_platform_counts <- expression_platform %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::group_by(which_platforms, platform) %>% 
  dplyr::summarise(total = sum(count)) %>% 
  dplyr::select(-platform) %>% 
  dplyr::distinct()

write.csv(expression_platform_counts, paste0(outdir, "/platform_boxplot_counts.csv"), quote = FALSE, row.names = FALSE)

# Stats
file.remove(paste0(outdir, "/platform_boxplot_stats.txt"))

for (platform_name in platform_labels) {
  
  expression_single_platform <- expression_platform %>% 
    dplyr::filter(platform == platform_name)
  
  kruskal.test(log2_cpm ~ which_platforms, expression_single_platform)
  wilcox_test <- pairwise.wilcox.test(expression_single_platform$log2_cpm, expression_single_platform$which_platforms, 
                                      p.adjust.method = "BH")
  
  sink(paste0(outdir, "/platform_boxplot_stats.txt"), append = TRUE)
  print(platform_name)
  print(wilcox_test)
  sink()
  
}

