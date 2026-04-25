#!/usr/bin/env Rscript
# H3K9me2 vs Gene Expression Correlation Analysis
# Window: 5000bp upstream, 500bp downstream of TSS
# Load required libraries
library(tidyverse)
library(ggplot2)
library(corrplot)
library(ggpubr)
# Set working directory and create output folder
setwd("/N/project/Krolab/isabella/ds-analysis")
dir.create("correlation_analysis_5kb", showWarnings = FALSE)

# 1. LOAD DATA
# Load peak-to-gene assignments for 5kb window
# This file has no headers, so we need to specify column names
peaks_near_tss <- read_tsv("peaks_near_TSS_5000_500.tsv", 
                           col_names = c("peak_chrom", "peak_start", "peak_end",
                                        "direction", "log2FC", "dot1",
                                        "gene_chrom", "gene_start", "gene_end", 
                                        "gene_name", "dot2", "strand")) %>%
  # Calculate distance to TSS (gene_start is the TSS for + strand genes)
  # use gene_start as TSS for all genes
  mutate(
    peak_midpoint = (peak_start + peak_end) / 2,
    distance_to_tss = peak_midpoint - gene_start
  ) %>%
  # Remove the "inf" values in log2FC
  filter(is.finite(log2FC))
# Load normalized gene expression data
gene_expression <- read_csv("logCPM_TMM_cleaned.csv") %>%
  # Strip version numbers from ensembl IDs
  mutate(ensembl_id = str_replace(GeneID, "\\.\\d+$", ""))
# Load gene ID mapping to convert between Ensembl IDs and gene symbols
gene_mapping <- read_tsv("/N/project/Krolab/isabella/annotations/ensembl_id_to_name.tsv",
                        col_names = c("ensembl_id_versioned", "gene_symbol")) %>%
  # Also strip version numbers from the mapping file
  mutate(ensembl_id = str_replace(ensembl_id_versioned, "\\.\\d+$", ""))
cat(sprintf("Loaded %d peak-gene associations\n", nrow(peaks_near_tss)))
cat(sprintf("Loaded %d genes with expression data\n", nrow(gene_expression)))

# 2. QUANTIFY H3K9me2 SIGNAL PER GENE
#  log2FC values are already in peaks_near_tss from bedtools intersect, no need to join with DiffReps output

# Calculate H3K9me2 metrics per gene
h3k9me2_per_gene <- peaks_near_tss %>%
  group_by(gene_name) %>%
  summarise(
    # Signal metrics using the log2FC from your file
    sum_log2FC = sum(log2FC, na.rm = TRUE),          # Total signal
    max_log2FC = max(log2FC, na.rm = TRUE),          # Strongest peak
    mean_log2FC = mean(log2FC, na.rm = TRUE),        # Average signal
    # Peak count
    n_peaks = n(),
    # Distance metrics
    min_distance = min(abs(distance_to_tss), na.rm = TRUE),
    mean_distance = mean(abs(distance_to_tss), na.rm = TRUE),
    # Track direction of regulation
    n_up = sum(direction == "Up"),
    n_down = sum(direction == "Down"),
    .groups = "drop"
  ) %>%
  # Replace infinite values with NA (only for double columns, not integers)
  mutate(across(where(is.double), ~na_if(., Inf))) %>%
  mutate(across(where(is.double), ~na_if(., -Inf)))

# 3. MERGE WITH GENE EXPRESSION DATA
# First, convert expression data to use gene symbols
# Join expression with mapping to get gene symbols
expression_with_symbols <- gene_expression %>%
  left_join(gene_mapping %>% select(ensembl_id, gene_symbol), 
            by = "ensembl_id") %>%
  filter(!is.na(gene_symbol))  # Remove genes without symbol mapping
cat(sprintf("Mapped %d genes from Ensembl IDs to gene symbols\n", nrow(expression_with_symbols)))
# Calculate mean expression per condition
expression_summary <- expression_with_symbols %>%
  rowwise() %>%
  mutate(
    control_mean = mean(c_across(starts_with("cntrl_rna")), na.rm = TRUE),
    SCZ_mean = mean(c_across(starts_with("scz_rna")), na.rm = TRUE),
    overall_mean = mean(c_across(starts_with("cntrl_rna") | starts_with("scz_rna")), na.rm = TRUE)
  ) %>%
  select(ensembl_id, gene_symbol, control_mean, SCZ_mean, overall_mean)
# Merge H3K9me2 signal with expression using gene symbols
merged_data <- h3k9me2_per_gene %>%
  inner_join(expression_summary, by = c("gene_name" = "gene_symbol")) %>%
  filter(!is.na(overall_mean))  # Remove genes with missing expression
cat(sprintf("Final dataset: %d genes with both H3K9me2 and expression data\n", nrow(merged_data)))
# Save merged dataset
write_csv(merged_data, "correlation_analysis_5kb/h3k9me2_expression_merged_5kb.csv")


# 4. CORRELATION ANALYSIS
# Function to calculate correlation with p-value
calc_correlation <- function(x, y, method = "pearson") {
  test <- cor.test(x, y, method = method)
  return(data.frame(
    correlation = test$estimate,
    p_value = test$p.value,
    method = method
  ))
}
# Correlation for different metrics and conditions
correlations <- list(
  # Using sum of log2FC
  sum_overall_pearson = calc_correlation(merged_data$sum_log2FC, 
                                         merged_data$overall_mean, "pearson"),
  sum_overall_spearman = calc_correlation(merged_data$sum_log2FC, 
                                          merged_data$overall_mean, "spearman"),
  sum_control_pearson = calc_correlation(merged_data$sum_log2FC, 
                                         merged_data$control_mean, "pearson"),
  sum_SCZ_pearson = calc_correlation(merged_data$sum_log2FC, 
                                     merged_data$SCZ_mean, "pearson"),
  # Using max log2FC
  max_overall_pearson = calc_correlation(merged_data$max_log2FC, 
                                         merged_data$overall_mean, "pearson"),
  max_overall_spearman = calc_correlation(merged_data$max_log2FC, 
                                          merged_data$overall_mean, "spearman"),
  # Using mean log2FC
  mean_overall_pearson = calc_correlation(merged_data$mean_log2FC, 
                                          merged_data$overall_mean, "pearson"),
  mean_overall_spearman = calc_correlation(merged_data$mean_log2FC, 
                                           merged_data$overall_mean, "spearman")
)
# Combine into a summary table
correlation_summary <- bind_rows(correlations, .id = "analysis") %>%
  mutate(
    significant = p_value < 0.05,
    FDR = p.adjust(p_value, method = "fdr")
  )
print(correlation_summary)
# Save correlation results
write_csv(correlation_summary, "correlation_analysis_5kb/correlation_summary_5kb.csv")


# 5. STRATIFIED ANALYSIS
# Bin genes by H3K9me2 signal strength
merged_data <- merged_data %>%
  mutate(
    h3k9me2_category = cut(sum_log2FC, 
                           breaks = quantile(sum_log2FC, probs = c(0, 0.25, 0.5, 0.75, 1)),
                           labels = c("Low", "Medium", "High", "Very High"),
                           include.lowest = TRUE)
  )

# 6. TOP CORRELATED GENES
# Calculate per-gene correlation might not be meaningful with limited samples
# Instead, identify genes with extreme H3K9me2 and expression patterns

# Top genes with high H3K9me2 and low expression (expected pattern)
high_h3k9me2_low_expr <- merged_data %>%
  filter(sum_log2FC > quantile(sum_log2FC, 0.75)) %>%
  filter(overall_mean < quantile(overall_mean, 0.25)) %>%
  arrange(desc(sum_log2FC)) %>%
  head(20)
write_csv(high_h3k9me2_low_expr, 
          "correlation_analysis_5kb/high_h3k9me2_low_expression_genes.csv")
# Top genes with low H3K9me2 and high expression
low_h3k9me2_high_expr <- merged_data %>%
  filter(sum_log2FC < quantile(sum_log2FC, 0.25)) %>%
  filter(overall_mean > quantile(overall_mean, 0.75)) %>%
  arrange(overall_mean) %>%
  head(20)
write_csv(low_h3k9me2_high_expr, 
          "correlation_analysis_5kb/low_h3k9me2_high_expression_genes.csv")

# 7. SUMMARY STATISTICS
summary_stats <- merged_data %>%
  summarise(
    n_genes = n(),
    mean_peaks_per_gene = mean(n_peaks),
    median_peaks_per_gene = median(n_peaks),
    mean_h3k9me2_signal = mean(sum_log2FC),
    sd_h3k9me2_signal = sd(sum_log2FC),
    mean_expression = mean(overall_mean),
    sd_expression = sd(overall_mean),
    min_distance_avg = mean(min_distance),
    max_distance_avg = mean(mean_distance)
  )

write_csv(summary_stats, "correlation_analysis_5kb/summary_statistics_5kb.csv")
print(summary_stats)
