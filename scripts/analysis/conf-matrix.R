#!/usr/bin/env Rscript

Sys.setenv(OPENBLAS_NUM_THREADS = 1)

library(data.table)
library(ggplot2)

cat("=== Creating Contingency Matrix with ALL DMGs ===\n")

# Load differential peaks
peaks <- fread("/N/project/Krolab/isabella/H3K9me2-Research/chip-seq/peak-analy/diffreps_output.txt")
cat("Total differential peaks:", nrow(peaks), "\n")

# Filter significant peaks (|log2FC| > 0.3) - KEEP Event column!
sig_peaks <- peaks[abs(log2FC) > 0.3, .(Chrom, Start, End, Event, log2FC)]
cat("Significant peaks (|log2FC| > 0.3):", nrow(sig_peaks), "\n")

# Load ALL peak-to-gene mappings
peak_gene <- fread("/N/project/Krolab/isabella/H3K9me2-Research/chip-seq/peak-analy/peak_to_gene.tsv",
                   header = FALSE)

# Assign column names
colnames(peak_gene) <- c("peak_chr", "peak_start", "peak_end", "peak_name", "peak_score", 
                         "strand1", "gene_chr", "gene_start", "gene_end", "gene_name", 
                         "dot", "strand2", "distance")

# Create peak keys
sig_peaks[, peak_key := paste(Chrom, Start, End, sep = "_")]
peak_gene[, peak_key := paste(peak_chr, peak_start, peak_end, sep = "_")]

# Merge significant peaks with gene names
dmg_data <- merge(sig_peaks, 
                  peak_gene[, .(peak_key, gene_name)], 
                  by = "peak_key", 
                  all.x = TRUE)

cat("Peak-gene associations:", nrow(dmg_data), "\n")

# Remove NAs
dmg_data <- dmg_data[!is.na(gene_name)]
cat("After removing NAs:", nrow(dmg_data), "\n")

# Get modification direction per gene (majority vote)
dmg_direction <- dmg_data[, .(
  scz_enriched_peaks = sum(Event == "Up"),
  control_enriched_peaks = sum(Event == "Down"),
  total_peaks = .N
), by = gene_name]

# Classify based on majority
dmg_direction[, h3k9me2_direction := ifelse(scz_enriched_peaks > control_enriched_peaks, 
                                             "SCZ-enriched", "Control-enriched")]

cat("\nUnique DMGs:", nrow(dmg_direction), "\n")
cat("DMGs by H3K9me2 direction:\n")
cat("  SCZ-enriched:", sum(dmg_direction$h3k9me2_direction == "SCZ-enriched"), "\n")
cat("  Control-enriched:", sum(dmg_direction$h3k9me2_direction == "Control-enriched"), "\n")

# Load expression data
expr <- fread("/N/project/Krolab/isabella/H3K9me2-Research/ds-analysis/logCPM_TMM_cleaned.csv")

control_cols <- grep("cntrl", colnames(expr), value = TRUE)
scz_cols <- grep("scz", colnames(expr), value = TRUE)

expr[, control_mean := rowMeans(.SD), .SDcols = control_cols]
expr[, scz_mean := rowMeans(.SD), .SDcols = scz_cols]
expr[, expression_change := scz_mean - control_mean]

# Load gene mapping
gene_map <- fread("/N/project/Krolab/isabella/H3K9me2-Research/annotations/ensembl_id_to_name.tsv",
                  header = FALSE, col.names = c("ensembl_id", "gene_name"))

# Merge expression with gene names
expr_with_names <- merge(expr, gene_map, by.x = "GeneID", by.y = "ensembl_id", all.x = TRUE)

# Merge DMG direction with expression
contingency_data <- merge(dmg_direction, 
                          expr_with_names[, .(gene_name, expression_change)],
                          by = "gene_name",
                          all.x = TRUE)

# Remove genes without expression data
contingency_data <- contingency_data[!is.na(expression_change)]

cat("\nDMGs with expression data:", nrow(contingency_data), "\n")

# Classify expression direction
contingency_data[, expression_direction := ifelse(expression_change > 0, 
                                                   "Higher in SCZ", "Higher in Control")]

cat("\nExpression direction distribution:\n")
cat("  Higher in SCZ:", sum(contingency_data$expression_direction == "Higher in SCZ"), "\n")
cat("  Higher in Control:", sum(contingency_data$expression_direction == "Higher in Control"), "\n")

# Create contingency table
contingency_table <- table(contingency_data$h3k9me2_direction, 
                           contingency_data$expression_direction)

cat("\n=== Contingency Table ===\n")
print(contingency_table)
cat("\nTotal genes:", sum(contingency_table), "\n")

# Create data frame for plotting
ct_df <- as.data.frame(contingency_table)
colnames(ct_df) <- c("H3K9me2", "Expression", "Count")

# Calculate percentages
ct_df$Percentage <- ct_df$Count / sum(ct_df$Count) * 100

# Chi-square test
chisq_result <- chisq.test(contingency_table)
cat("\nChi-square test:\n")
cat("  X-squared:", chisq_result$statistic, "\n")
cat("  p-value:", chisq_result$p.value, "\n")

# Calculate consistency
if (nrow(contingency_table) >= 2 && ncol(contingency_table) >= 2) {
  consistent <- contingency_table["Control-enriched", "Higher in Control"] + 
                contingency_table["SCZ-enriched", "Higher in SCZ"]
  consistency_pct <- consistent / sum(contingency_table) * 100
} else {
  consistency_pct <- NA
}

cat("Consistency with repression:", round(consistency_pct, 1), "%\n")

# Reorder for plotting
ct_df$H3K9me2 <- factor(ct_df$H3K9me2, levels = c("SCZ-enriched", "Control-enriched"))
ct_df$Expression <- factor(ct_df$Expression, levels = c("Higher in Control", "Higher in SCZ"))

# Create heatmap
p <- ggplot(ct_df, aes(x = Expression, y = H3K9me2, fill = Count)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = paste0(Count, "\n(", round(Percentage, 1), "%)")), 
            color = "white", size = 8, fontface = "bold") +
  scale_fill_gradient(low = "#5F9EA0", high = "#2F4F4F", name = "# genes") +
  labs(title = "H3K9me2 Modification vs Expression Direction",
       subtitle = paste0("Chi-square p = ", round(chisq_result$p.value, 3), 
                        " | Consistent with repression: ", round(consistency_pct, 1), "%"),
       x = "Differential Expression",
       y = "Differential H3K9me2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave("contingency_matrix_ALL_DMGs.png", p, width = 10, height = 6, dpi = 300, bg = "white")

cat("\n=== DONE! ===\n")
cat("Saved: contingency_matrix_ALL_DMGs.png\n")
cat("\nFinal counts:\n")
print(contingency_table)
