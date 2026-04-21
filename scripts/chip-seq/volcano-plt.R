# Load data
df <- read.table("/N/project/Krolab/isabella/chip-seq/peak-analy/diffreps_output.txt",
                 header=TRUE, sep="\t", comment.char="#")

# Calculate -log10(p-value)
df$log10p <- -log10(df$pval)

# Classify significance
df$EventColor <- ifelse(df$padj < 0.05 & abs(df$log2FC) > 1,
                        ifelse(df$log2FC > 1, "Up", "Down"), "NS")

# Count for annotation
num_up <- sum(df$EventColor == "Up", na.rm=TRUE)
num_down <- sum(df$EventColor == "Down", na.rm=TRUE)
num_sig <- num_up + num_down

# Save plot
png("/N/project/Krolab/isabella/chip-seq/peak-analy/volcano_plot.png", width=1000, height=800)
par(mar = c(5, 5, 4, 6))  # Set margins to avoid text cutoff

# Main plot
plot(df$log2FC, df$log10p,
     col=ifelse(df$EventColor == "Up", "red",
                ifelse(df$EventColor == "Down", "blue", "gray")),
     pch=20,
     main="Volcano Plot of Differential H3K9me2 Peaks",
     xlab="Log2 Fold Change", ylab="-log10(p-value)",
     cex.main=2, cex.lab=1.6, cex.axis=1.4)

# Add threshold lines
abline(h = -log10(0.05), col="black", lty=2, lwd=1.5)  # Horizontal line at p=0.05
abline(v = -2, col="black", lty=2, lwd=1.5)            # Vertical lines at -2 and 2
abline(v = 2, col="black", lty=2, lwd=1.5)

# Add counts as text on plot
legend("topright", legend=c("Up", "Down", "NS"),
       col=c("red", "blue", "gray"), pch=20, cex=1.4)

text(x = min(df$log2FC, na.rm=TRUE), y = max(df$log10p, na.rm=TRUE),
     labels = paste0("Up: ", num_up, "\nDown: ", num_down, "\nTotal sig: ", num_sig),
     adj = c(0,1), cex=1.4)

dev.off()

