#!/usr/bin/env Rscript
# ============================================================================
# Per-Sample Correlation Boxplot: Control vs SCZ
# Calculates correlation between H3K9me2 and expression for EACH sample
# ============================================================================

library(tidyverse)
library(ggplot2)

setwd("/N/project/Krolab/isabella/ds-analysis")
dir.create("correlation_analysis_5kb", showWarnings = FALSE)

cat("Loading merged data...\n")

# ============================================================================
# 1. LOAD THE MERGED DATA
# ============================================================================

# Load the merged dataset we created earlier
merged_data <- read_csv("correlation_analysis_5kb/h3k9me2_expression_merged_5kb.csv")

# Load the full expression data (we need individual sample columns)
gene_expression <- read_csv("logCPM_TMM_cleaned.csv") %>%
  mutate(ensembl_id = str_replace(GeneID, "\\.\\d+$", ""))

# Load gene mapping
gene_mapping <- read_tsv("/N/project/Krolab/isabella/annotations/ensembl_id_to_name.tsv",
                        col_names = c("ensembl_id_versioned", "gene_symbol")) %>%
  mutate(ensembl_id = str_replace(ensembl_id_versioned, "\\.\\d+$", ""))

# Add gene symbols to expression data
expression_with_symbols <- gene_expression %>%
  left_join(gene_mapping %>% select(ensembl_id, gene_symbol), by = "ensembl_id") %>%
  filter(!is.na(gene_symbol))

# ============================================================================
# 2. PREPARE DATA FOR PER-SAMPLE CORRELATIONS
# ============================================================================

cat("Preparing data for per-sample correlations...\n")

# Keep only genes that have H3K9me2 peaks (those in merged_data)
genes_with_peaks <- unique(merged_data$gene_name)

# Filter expression data to only these genes
expression_filtered <- expression_with_symbols %>%
  filter(gene_symbol %in% genes_with_peaks) %>%
  select(gene_symbol, starts_with("cntrl_rna"), starts_with("scz_rna"))

# Get H3K9me2 signal per gene (we'll use sum_log2FC)
h3k9me2_signal <- merged_data %>%
  select(gene_name, sum_log2FC)

# Merge
data_for_correlation <- expression_filtered %>%
  left_join(h3k9me2_signal, by = c("gene_symbol" = "gene_name"))

cat(sprintf("Analyzing %d genes across all samples\n", nrow(data_for_correlation)))

# ============================================================================
# 3. CALCULATE CORRELATION FOR EACH SAMPLE
# ============================================================================

cat("\nCalculating per-sample correlations...\n")

# Get all sample column names
control_samples <- colnames(data_for_correlation)[str_detect(colnames(data_for_correlation), "^cntrl_rna")]
scz_samples <- colnames(data_for_correlation)[str_detect(colnames(data_for_correlation), "^scz_rna")]

# Function to calculate correlation for one sample
calc_sample_correlation <- function(sample_col, data, method = "pearson") {
  x <- data$sum_log2FC
  y <- data[[sample_col]]
  
  # Remove NA values
  valid_idx <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
  
  if (sum(valid_idx) < 3) {
    return(NA)
  }
  
  cor(x[valid_idx], y[valid_idx], method = method)
}

# Calculate correlations for all control samples
control_correlations <- tibble(
  sample = control_samples,
  correlation = -map_dbl(control_samples, ~calc_sample_correlation(.x, data_for_correlation)),  # NEGATIVE
  group = "Control"
)

# Calculate correlations for all SCZ samples
scz_correlations <- tibble(
  sample = scz_samples,
  correlation = -map_dbl(scz_samples, ~calc_sample_correlation(.x, data_for_correlation)),  # NEGATIVE
  group = "Case"  # Changed from "SCZ" to "Case" to match your plot
)

# Combine
all_correlations <- bind_rows(control_correlations, scz_correlations) %>%
  filter(!is.na(correlation)) %>%
  # Set factor levels to control order (Control on left, Case on right)
  mutate(group = factor(group, levels = c("Control", "Case")))

cat(sprintf("Control samples: %d correlations calculated\n", sum(all_correlations$group == "Control")))
cat(sprintf("Case samples: %d correlations calculated\n", sum(all_correlations$group == "Case")))

# Print summary statistics
cat("\nSummary statistics:\n")
summary_stats <- all_correlations %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(correlation),
    median = median(correlation),
    sd = sd(correlation),
    min = min(correlation),
    max = max(correlation)
  )
print(summary_stats)

# Save the correlation data
write_csv(all_correlations, "correlation_analysis_5kb/per_sample_correlations.csv")

# ============================================================================
# 4. CREATE THE BOXPLOT
# ============================================================================

cat("\nCreating boxplot...\n")

# Create the plot matching the style of the uploaded image
p <- ggplot(all_correlations, aes(x = group, y = correlation, fill = group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Hide outliers from boxplot
  geom_jitter(width = 0.2, alpha = 0.5, size = 2.5, color = "gray40") +  # Add individual points
  scale_fill_manual(values = c("Control" = "#7AB8A8", "Case" = "#E89B7E")) +  # Match colors
  labs(
    title = "Pearson Correlation by Group",
    x = "Sample Group",
    y = "Negative Pearson Correlation (-r)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )

ggsave("correlation_analysis_5kb/boxplot_per_sample_correlations.png", 
       p, width = 8, height = 6, dpi = 300)

cat("\nBoxplot saved to: correlation_analysis_5kb/boxplot_per_sample_correlations.png\n")

# ============================================================================
# 5. STATISTICAL TEST
# ============================================================================

cat("\nPerforming statistical test...\n")

# Test if correlations differ between groups
control_cors <- all_correlations %>% filter(group == "Control") %>% pull(correlation)
case_cors <- all_correlations %>% filter(group == "Case") %>% pull(correlation)

# Wilcoxon test (non-parametric)
wilcox_result <- wilcox.test(control_cors, case_cors)
cat(sprintf("Wilcoxon rank-sum test p-value: %.4f\n", wilcox_result$p.value))

# T-test (parametric)
t_result <- t.test(control_cors, case_cors)
cat(sprintf("T-test p-value: %.4f\n", t_result$p.value))

# Save test results
test_results <- tibble(
  test = c("Wilcoxon", "T-test"),
  p_value = c(wilcox_result$p.value, t_result$p.value),
  control_mean = c(mean(control_cors), mean(control_cors)),
  case_mean = c(mean(case_cors), mean(case_cors)),
  difference = c(mean(control_cors) - mean(case_cors), mean(control_cors) - mean(case_cors))
)

write_csv(test_results, "correlation_analysis_5kb/group_comparison_stats.csv")

# ============================================================================
# COMPLETE
# ============================================================================

cat("\n=== Analysis Complete ===\n")
cat("Generated files:\n")
cat("  - boxplot_per_sample_correlations.png (main plot)\n")
cat("  - per_sample_correlations.csv (correlation values)\n")
cat("  - group_comparison_stats.csv (statistical test results)\n")
