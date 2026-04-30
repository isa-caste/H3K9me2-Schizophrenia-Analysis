#!/usr/bin/env Rscript
# Gene Ontology Enrichment Analysis
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
.libPaths(c("/N/project/Krolab/isabella/H3K9me2-Research/R_libs", .libPaths()))
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
setwd("/N/project/Krolab/isabella/H3K9me2-Research/ds-analysis")
cat("=== Gene Ontology Enrichment Analysis ===\n")
cat("Filtering for: H3K9me2↑ & Expression↓ coordinated changes\n\n")
# 0. USER-DEFINED THRESHOLDS
cat("0. Setting filtering thresholds...\n")
EXPR_LOGFC_THRESHOLD <- -0.5
EXPR_PADJ_THRESHOLD <- 0.05
GO_PADJ_THRESHOLD <- 0.05
GO_QVALUE_THRESHOLD <- 0.2
cat(sprintf("  Expression logFC < %.2f, padj < %.3f\n", EXPR_LOGFC_THRESHOLD, EXPR_PADJ_THRESHOLD))
cat(sprintf("  GO/KEGG enrichment padj < %.3f\n\n", GO_PADJ_THRESHOLD))
# 1. LOAD DATA
cat("1. Loading data...\n")
raw_data <- read.csv("peak_gene_DEGs_unfiltered.tsv", sep = "\t", stringsAsFactors = FALSE)
cat(sprintf("  Loaded: %d rows, %d columns\n", nrow(raw_data), ncol(raw_data)))
# 2. PREPARE DATA
cat("\n2. Preparing data for filtering...\n")
data <- data.frame(
  Gene_Name = raw_data$gene_name,
  Ensembl_ID = raw_data$gene_id,
  expr_logFC = raw_data$logFC,
  expr_padj = raw_data$adj.P.Val,
  expr_praw = raw_data$P.Value,
  expr_downregulated = raw_data$is_downregulated,
  peak_direction = raw_data$peak_direction,
  peak_significant = raw_data$is_significant_adj,
  stringsAsFactors = FALSE
)
data <- data[!is.na(data$Gene_Name) & !is.na(data$expr_logFC), ]
data$peak_significant <- data$peak_significant == "True"
data$expr_downregulated <- data$expr_downregulated == "True"
cat(sprintf("  Prepared: %d genes with complete data\n", nrow(data)))
cat("\n  Summary of expression changes:\n")
cat(sprintf("    Downregulated genes: %d\n", sum(data$expr_downregulated == TRUE, na.rm = TRUE)))
cat(sprintf("    Peak direction Up: %d\n", sum(data$peak_direction == "Up", na.rm = TRUE)))
cat(sprintf("    Peak direction Down: %d\n", sum(data$peak_direction == "Down", na.rm = TRUE)))
cat(sprintf("    Peak significant TRUE: %d\n", sum(data$peak_significant == TRUE, na.rm = TRUE)))
# 3. FILTER FOR COORDINATED CHANGES
cat("\n3. Filtering for coordinated regulatory changes...\n")
filter_idx <- (data$peak_direction == "Up") &
              (data$expr_logFC < EXPR_LOGFC_THRESHOLD) &
	      (data$expr_praw < EXPR_PADJ_THRESHOLD) &
              (data$expr_downregulated == TRUE)
coordinated_changes <- data[filter_idx, ]
coordinated_changes <- coordinated_changes[!duplicated(coordinated_changes$Gene_Name), ]
cat(sprintf("  Genes with coordinated silencing (H3K9me2↑ & Expression↓): %d\n", nrow(coordinated_changes)))
if (nrow(coordinated_changes) > 0) {
  cat("\n  Top 10 filtered genes:\n")
  print(head(coordinated_changes[, c("Gene_Name", "expr_logFC", "expr_padj", "peak_direction")], 10))
}
all_genes_with_peaks <- unique(data$Gene_Name[data$peak_significant == TRUE])
cat(sprintf("\n  Background gene lists:\n"))
# 5. SAVE FILTERED GENE LIST
cat("\n5. Saving filtered gene lists...\n")
dir.create("enrichment_results", showWarnings = FALSE)
write.csv(coordinated_changes, "enrichment_results/genes_coordinated_changes_filtered.csv", row.names = FALSE)
cat(sprintf("  - genes_coordinated_changes_filtered.csv (%d genes)\n", nrow(coordinated_changes)))
writeLines(coordinated_changes$Gene_Name, "enrichment_results/gene_list_coordinated.txt")
cat(sprintf("  - gene_list_coordinated.txt\n"))
writeLines(all_genes_with_peaks, "enrichment_results/gene_list_background_peaks.txt")
cat(sprintf("  - gene_list_background_peaks.txt\n\n"))
cat(sprintf("    All genes with significant peaks: %d genes\n", length(all_genes_with_peaks)))

# 4. COMPARE WITH REFERENCE PAPER
cat("\n4. Gene overlap analysis with reference paper...\n")
PAPER_GENE_FILE <- "reference_paper_genes.txt"
if (file.exists(PAPER_GENE_FILE)) {
  paper_genes <- tolower(trimws(readLines(PAPER_GENE_FILE)))
  coordinated_genes_lower <- tolower(coordinated_changes$Gene_Name)
  confirmed_idx <- coordinated_genes_lower %in% paper_genes
  novel_idx <- !(coordinated_genes_lower %in% paper_genes)
  confirmed_genes <- coordinated_changes[confirmed_idx, ]
  if (nrow(confirmed_genes) > 0) confirmed_genes$Status <- "Confirmed (in paper)"
  novel_genes <- coordinated_changes[novel_idx, ]
  if (nrow(novel_genes) > 0) novel_genes$Status <- "Novel (not in paper)"
  cat(sprintf("  Confirmed genes (in paper): %d\n", nrow(confirmed_genes)))
  cat(sprintf("  Novel genes (not in paper): %d\n", nrow(novel_genes)))
  write.csv(confirmed_genes, "genes_confirmed_in_paper.csv", row.names = FALSE)
  write.csv(novel_genes, "genes_novel_findings.csv", row.names = FALSE)
} else {
  cat(sprintf("  WARNING: %s not found. Skipping overlap analysis.\n", PAPER_GENE_FILE))
  confirmed_genes <- coordinated_changes
  if (nrow(confirmed_genes) > 0) confirmed_genes$Status <- "Analyzed gene"
}
# 5. SAVE FILTERED GENE LIST
cat("\n5. Saving filtered gene lists...\n")
dir.create("enrichment_results", showWarnings = FALSE)
write.csv(coordinated_changes, "enrichment_results/genes_coordinated_changes_filtered.csv", row.names = FALSE)
cat(sprintf("  - genes_coordinated_changes_filtered.csv (%d genes)\n", nrow(coordinated_changes)))
writeLines(coordinated_changes$Gene_Name, "enrichment_results/gene_list_coordinated.txt")
cat(sprintf("  - gene_list_coordinated.txt\n"))
writeLines(all_genes_with_peaks, "enrichment_results/gene_list_background_peaks.txt")
cat(sprintf("  - gene_list_background_peaks.txt\n\n"))
# 6. CONVERT GENE SYMBOLS TO ENTREZ IDs
cat("6. Converting gene symbols to Entrez IDs...\n")
if (nrow(coordinated_changes) == 0) {
  cat("  ERROR: No genes passed filters. Cannot proceed with GO analysis.\n")
  cat("  Check that is_significant_adj contains 'True' values in your input data.\n")
  quit(save = "no", status = 1)
}
gene_list_entrez <- bitr(coordinated_changes$Gene_Name,
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db,
                         drop = TRUE)
background_entrez <- bitr(all_genes_with_peaks,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db,
                          drop = TRUE)
cat(sprintf("  Converted %d/%d genes to Entrez IDs\n", nrow(gene_list_entrez), nrow(coordinated_changes)))
cat(sprintf("  Background: %d/%d genes\n\n", nrow(background_entrez), length(all_genes_with_peaks)))
# 7. GO ENRICHMENT ANALYSIS
cat("7. Running GO enrichment analysis...\n")
cat("  - Biological Process...\n")
go_bp <- enrichGO(gene = gene_list_entrez$ENTREZID,
                  universe = background_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = GO_PADJ_THRESHOLD,
                  qvalueCutoff = GO_QVALUE_THRESHOLD,
                  readable = TRUE)
cat("  - Molecular Function...\n")
go_mf <- enrichGO(gene = gene_list_entrez$ENTREZID,
                  universe = background_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = GO_PADJ_THRESHOLD,
                  qvalueCutoff = GO_QVALUE_THRESHOLD,
                  readable = TRUE)
cat("  - Cellular Component...\n")
go_cc <- enrichGO(gene = gene_list_entrez$ENTREZID,
                  universe = background_entrez$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = GO_PADJ_THRESHOLD,
                  qvalueCutoff = GO_QVALUE_THRESHOLD,
                  readable = TRUE)
cat("  - Done!\n\n")
# 8. KEGG PATHWAY ANALYSIS
cat("8. Running KEGG pathway analysis...\n")
kegg <- enrichKEGG(gene = gene_list_entrez$ENTREZID,
                   universe = background_entrez$ENTREZID,
                   organism = "hsa",
                   pvalueCutoff = GO_PADJ_THRESHOLD,
                   qvalueCutoff = GO_QVALUE_THRESHOLD)
if (!is.null(kegg) && nrow(kegg) > 0){
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}
cat("  - Done!\n\n")
# 9. SAVE ENRICHMENT RESULTS
cat("9. Saving enrichment results...\n")
if (!is.null(go_bp) && nrow(go_bp) > 0) {
  write.csv(as.data.frame(go_bp), "enrichment_results/GO_Biological_Process.csv", row.names = FALSE)
  cat(sprintf("  - Found %d enriched BP terms\n", nrow(go_bp)))
} else {
  cat("  - No significant BP terms\n")
}
if (!is.null(go_mf) && nrow(go_mf) > 0) {
  write.csv(as.data.frame(go_mf), "enrichment_results/GO_Molecular_Function.csv", row.names = FALSE)
  cat(sprintf("  - Found %d enriched MF terms\n", nrow(go_mf)))
} else {
  cat("  - No significant MF terms\n")
}
if (!is.null(go_cc) && nrow(go_cc) > 0) {
  write.csv(as.data.frame(go_cc), "enrichment_results/GO_Cellular_Component.csv", row.names = FALSE)
  cat(sprintf("  - Found %d enriched CC terms\n", nrow(go_cc)))
} else {
  cat("  - No significant CC terms\n")
}
if (!is.null(kegg) && nrow(go_kegg) > 0) {
  write.csv(as.data.frame(kegg), "enrichment_results/KEGG_Pathways.csv", row.names = FALSE)
  cat(sprintf("  - Found %d enriched pathways\n\n", nrow(kegg)))
} else {
  cat("  - No significant pathways\n\n")
}
# 10. VISUALIZATIONS
cat("10. Creating visualizations...\n")
if (!is.null(go_bp) && nrow(go_bp) > 0) {
  p1 <- barplot(go_bp, showCategory = 20) +
    ggtitle("Top 20 Enriched Biological Processes\n(Genes with H3K9me2↑ & Expression↓)") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("enrichment_results/01_GO_BP_barplot.png", p1, width = 12, height = 10, dpi = 300)
  cat("  - GO BP barplot\n")
  p2 <- dotplot(go_bp, showCategory = 20) +
    ggtitle("GO Biological Process Enrichment") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("enrichment_results/02_GO_BP_dotplot.png", p2, width = 12, height = 10, dpi = 300)
  cat("  - GO BP dotplot\n")
}
if (!is.null(go_mf) && nrow(go_mf) > 0) {
  p3 <- barplot(go_mf, showCategory = 20) +
    ggtitle("Top 20 Enriched Molecular Functions") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("enrichment_results/03_GO_MF_barplot.png", p3, width = 12, height = 10, dpi = 300)
  cat("  - GO MF barplot\n")
  p3b <- dotplot(go_mf, showCategory = 20) +
    ggtitle("GO Molecular Function Enrichment") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("enrichment_results/03b_GO_MF_dotplot.png", p3b, width = 12, height = 10, dpi = 300)
}
if (!is.null(go_cc) && nrow(go_cc) > 0)  {
  p3c <- barplot(go_cc, showCategory = 15) +
    ggtitle("Top 15 Enriched Cellular Components") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("enrichment_results/04_GO_CC_barplot.png", p3c, width = 12, height = 8, dpi = 300)
  cat("  - GO CC barplot\n")
}
if (!is.null(kegg) && nrow(kegg) > 0){
  p4 <- dotplot(kegg, showCategory = 20) +
    ggtitle("KEGG Pathway Enrichment") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("enrichment_results/05_KEGG_dotplot.png", p4, width = 12, height = 10, dpi = 300)
  cat("  - KEGG dotplot\n")
  p4b <- barplot(kegg, showCategory = 20) +
    ggtitle("Top 20 KEGG Pathways") +
    theme(axis.text.y = element_text(size = 9))
  ggsave("enrichment_results/05b_KEGG_barplot.png", p4b, width = 12, height = 10, dpi = 300)
  cat("  - KEGG barplot\n")
}
cat("  - Done!\n\n")
# 11. SUMMARY REPORT
cat("11. Generating summary report...\n\n")
summary_data <- data.frame(
  Analysis = c("GO Biological Process", "GO Molecular Function", "GO Cellular Component", "KEGG Pathways"),
  Significant_Terms = c(nrow(go_bp), nrow(go_mf), nrow(go_cc), nrow(kegg)),
  Input_Genes = rep(nrow(gene_list_entrez), 4),
  Background_Genes = rep(nrow(background_entrez), 4),
  stringsAsFactors = FALSE
)
print(summary_data)
write.csv(summary_data, "enrichment_results/summary_statistics.csv", row.names = FALSE)
if (!is.null(go_bp) && nrow(go_bp) > 0) {
  cat("\n=== Top 10 Biological Processes ===\n")
  top_bp <- as.data.frame(go_bp)[1:min(10, nrow(go_bp)), c("Description", "GeneRatio", "pvalue", "p.adjust", "Count")]
  print(top_bp)
  write.csv(top_bp, "enrichment_results/top_BP_terms.csv", row.names = FALSE)
}
if (!is.null(go_mf) && nrow(go_mf) > 0) {
  cat("\n=== Top 10 Molecular Functions ===\n")
  top_mf <- as.data.frame(go_mf)[1:min(10, nrow(go_mf)), c("Description", "GeneRatio", "pvalue", "p.adjust", "Count")]
  print(top_mf)
  write.csv(top_mf, "enrichment_results/top_MF_terms.csv", row.names = FALSE)
}
if (!is.null(kegg) && nrow(kegg) > 0) {
  cat("\n=== Top 10 KEGG Pathways ===\n")
  top_kegg <- as.data.frame(kegg)[1:min(10, nrow(kegg)), c("Description", "GeneRatio", "pvalue", "p.adjust", "Count")]
  print(top_kegg)
  write.csv(top_kegg, "enrichment_results/top_KEGG_pathways.csv", row.names = FALSE)
}

