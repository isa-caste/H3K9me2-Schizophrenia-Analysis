#!/usr/bin/env Rscript

library(edgeR)

# Define input and output paths
input_path <- "/N/project/Krolab/isabella/rna-seq/quant-norm/featurecounts/featureCounts_output.txt"
output_path <- "/N/project/Krolab/isabella/rna-seq/quant-norm/logCPM_TMM.csv"

#  counts file
counts <- read.delim(input_path, comment.char = "#")

# Keep only geneid and sample count columns (remove Chr, Start, etc.)
counts_clean <- counts[ , c(1, 7:ncol(counts)) ]
rownames(counts_clean) <- counts_clean$Geneid
counts_clean <- counts_clean[ , -1 ]  # remove Geneid column

# create DGEList object
dge <- DGEList(counts = counts_clean)

# filter low-expression genes: CPM > 1 in at least 1 sample
keep <- rowSums(cpm(dge) > 1) >= 1
dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# log2 CPM with prior.count = 3
logCPM <- cpm(dge, log = TRUE, prior.count = 3)

# output
write.csv(logCPM, file = output_path)

