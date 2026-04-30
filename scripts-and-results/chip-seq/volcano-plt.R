# Load data
df <- read.table("/N/project/Krolab/isabella/H3K9me2-Research/chip-seq/peak-analy/diffreps_output.txt",
                 header=TRUE, sep="\t", comment.char="#")

# Calculate -log10(p-value)
df$log10p <- -log10(df$pval)

# Significance thresholds
fc_threshold <- 0.3     # log2 fold change cutoff
padj_threshold <- 0.05  # adjusted p-value cutoff

# Classify significance
df$EventColor <- ifelse(df$padj < padj_threshold & df$log2FC > fc_threshold, "Up",
                  ifelse(df$padj < padj_threshold & df$log2FC < -fc_threshold, "Down", "NS"))

# Count for annotation
num_up <- sum(df$EventColor == "Up", na.rm=TRUE)
num_down <- sum(df$EventColor == "Down", na.rm=TRUE)
num_sig <- num_up + num_down

# Save plot as PNG
png("/N/project/Krolab/isabella/H3K9me2-Research/chip-seq/peak-analy/volcano_plot.png",
    width=1000, height=800, type="cairo")

par(mar = c(5, 5, 4, 6))

# Main plot
plot(df$log2FC, df$log10p,
     col=ifelse(df$EventColor == "Up", "red",
                ifelse(df$EventColor == "Down", "blue", "gray")),
     pch=20,
     main="Volcano Plot of Differential H3K9me2 Peaks",
     xlab="Log2 Fold Change", ylab="-log10(p-value)",
     cex.main=2, cex.lab=1.6, cex.axis=1.4)

# Add threshold lines
abline(h = -log10(padj_threshold), col="black", lty=2, lwd=1.5)
abline(v = -fc_threshold, col="black", lty=2, lwd=1.5)
abline(v = fc_threshold, col="black", lty=2, lwd=1.5)

# Legend
legend("topright", legend=c("Up", "Down", "NS"),
       col=c("red", "blue", "gray"), pch=20, cex=1.4)

# Annotation with counts
text(x = min(df$log2FC, na.rm=TRUE), y = max(df$log10p, na.rm=TRUE),
     labels = paste0("Up: ", num_up, "\nDown: ", num_down, "\nTotal sig: ", num_sig),
     adj = c(0,1), cex=1.4)

dev.off()
