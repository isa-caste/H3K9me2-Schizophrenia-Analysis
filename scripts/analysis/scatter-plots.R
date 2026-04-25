#!/usr/bin/env Rscript
# H3K9me2 Peak COUNT vs Gene Expression Scatter Plots
# Using number of peaks instead of signal strength


library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd("/N/project/Krolab/isabella/ds-analysis")

cat("Loading merged data...\n")

# Load the merged dataset (already has n_peaks calculated)
merged_data <- read_csv("correlation_analysis_5kb/h3k9me2_expression_merged_5kb.csv")

cat(sprintf("Loaded %d genes\n", nrow(merged_data)))
cat(sprintf("Peak count range: %d to %d\n", min(merged_data$n_peaks), max(merged_data$n_peaks)))
cat(sprintf("Mean peaks per gene: %.1f\n", mean(merged_data$n_peaks)))
cat(sprintf("Median peaks per gene: %.0f\n", median(merged_data$n_peaks)))


# OPTIONAL FILTERING
threshold_min_peaks <- 0  # Keep genes with at least this many peaks
                        
# use log scale for better visualization
use_log_scale <- TRUE 
filtered_data <- merged_data %>%
  filter(n_peaks >= threshold_min_peaks)
cat(sprintf("\nFiltered to genes with >= %d peaks: %d genes (%.1f%% of original)\n",
            threshold_min_peaks, nrow(filtered_data), 
            100 * nrow(filtered_data) / nrow(merged_data)))


# CALCULATE CORRELATIONS
cat("\nCalculating correlations...\n")
cor_overall <- cor.test(filtered_data$n_peaks, filtered_data$overall_mean, method = "pearson")
cat(sprintf("Peak count vs Overall expression: R = %.3f, p = %.2e\n", 
            cor_overall$estimate, cor_overall$p.value))
cor_control <- cor.test(filtered_data$n_peaks, filtered_data$control_mean, method = "pearson")
cat(sprintf("Peak count vs Control expression: R = %.3f, p = %.2e\n", 
            cor_control$estimate, cor_control$p.value))
cor_scz <- cor.test(filtered_data$n_peaks, filtered_data$SCZ_mean, method = "pearson")
cat(sprintf("Peak count vs SCZ expression: R = %.3f, p = %.2e\n", 
            cor_scz$estimate, cor_scz$p.value))


# PLOT 1: Peak Count vs Overall Expression
p1 <- ggplot(filtered_data, aes(x = n_peaks, y = overall_mean)) +
  geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = 0.7, label.y.npc = 0.9) +
  labs(
    title = "H3K9me2 Peak Count vs Gene Expression (5kb window)",
    subtitle = sprintf("n = %d genes (>= %d peaks)", nrow(filtered_data), threshold_min_peaks),
    x = "Number of H3K9me2 Peaks (5kb window)",
    y = "Mean Gene Expression (log2 CPM)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

# Add log scale
if (use_log_scale) {
  p1 <- p1 + scale_x_log10(breaks = c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000))
}
ggsave("correlation_analysis_5kb/scatter_peak_count_vs_expression.png",
       p1, width = 8, height = 6, dpi = 300)


# PLOT 2: Peak Count vs Expression - Control vs SCZ Comparison
p2 <- ggplot(filtered_data, aes(x = n_peaks, y = control_mean)) +
  geom_point(alpha = 0.5, size = 2, color = "blue") +
  geom_smooth(method = "lm", color = "darkblue", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = 0.7, label.y.npc = 0.9) +
  labs(
    title = "Control: H3K9me2 Peak Count vs Expression",
    x = "Number of H3K9me2 Peaks",
    y = "Mean Expression (log2 CPM)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

p3 <- ggplot(filtered_data, aes(x = n_peaks, y = SCZ_mean)) +
  geom_point(alpha = 0.5, size = 2, color = "red") +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  stat_cor(method = "pearson", label.x.npc = 0.7, label.y.npc = 0.9) +
  labs(
    title = "SCZ: H3K9me2 Peak Count vs Expression",
    x = "Number of H3K9me2 Peaks",
    y = "Mean Expression (log2 CPM)"
  ) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold"))

# Add log scale
if (use_log_scale) {
  p2 <- p2 + scale_x_log10(breaks = c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000))
  p3 <- p3 + scale_x_log10(breaks = c(1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000))
}
p_combined <- ggarrange(p2, p3, ncol = 2, nrow = 1)
ggsave("correlation_analysis_5kb/scatter_peak_count_control_vs_SCZ.png",
       p_combined, width = 12, height = 5, dpi = 300)

# SUMMARY STATISTICS
cat("\nGenerating summary statistics...\n")
# Peak count bins
peak_count_summary <- filtered_data %>%
  mutate(
    peak_bin = case_when(
      n_peaks == 1 ~ "1 peak",
      n_peaks == 2 ~ "2 peaks",
      n_peaks == 3 ~ "3 peaks",
      n_peaks >= 4 & n_peaks <= 5 ~ "4-5 peaks",
      n_peaks >= 6 & n_peaks <= 10 ~ "6-10 peaks",
      n_peaks > 10 ~ ">10 peaks"
    )
  ) %>%
  group_by(peak_bin) %>%
  summarise(
    n_genes = n(),
    mean_expression = mean(overall_mean),
    median_expression = median(overall_mean),
    .groups = "drop"
  ) %>%
  arrange(factor(peak_bin, levels = c("1 peak", "2 peaks", "3 peaks", "4-5 peaks", "6-10 peaks", ">10 peaks")))
print(peak_count_summary)
write_csv(peak_count_summary, "correlation_analysis_5kb/peak_count_summary.csv")
# Overall summary
overall_summary <- tibble(
  metric = c("Total genes", "Min peaks", "Max peaks", "Mean peaks", "Median peaks",
             "Correlation (Pearson)", "P-value"),
  value = c(
    as.character(nrow(filtered_data)),
    as.character(min(filtered_data$n_peaks)),
    as.character(max(filtered_data$n_peaks)),
    sprintf("%.2f", mean(filtered_data$n_peaks)),
    as.character(median(filtered_data$n_peaks)),
    sprintf("%.3f", cor_overall$estimate),
    sprintf("%.2e", cor_overall$p.value)
  )
)
write_csv(overall_summary, "correlation_analysis_5kb/peak_count_overall_summary.csv")
print(overall_summary)

cat("\n=== Analysis Complete ===\n")
cat("Generated files:\n")
cat("  - scatter_peak_count_vs_expression.png (main plot)\n")
cat("  - scatter_peak_count_control_vs_SCZ.png (comparison)\n")
cat("  - peak_count_summary.csv (stats by peak count bins)\n")
cat("  - peak_count_overall_summary.csv (overall statistics)\n")
cat(sprintf("\nUsing genes with >= %d peaks\n", threshold_min_peaks))
