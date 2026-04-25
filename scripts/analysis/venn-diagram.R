# DMG vs DEG Venn Diagram 
library(ggplot2)

PEAK_GENE_FILE <- "peak_gene_DEGs_unfiltered.tsv"
LIMMA_FILE     <- "/N/project/Krolab/isabella/H3K9me2-Research/rna-seq/DEG/limma_deg_results.csv"
COORD_FILE     <- "enrichment_results/genes_coordinated_changes_filtered.csv"
ENSEMBL_MAP    <- "/N/project/Krolab/isabella/H3K9me2-Research/annotations/ensembl_id_to_name.tsv"
OUT_DIR        <- "/N/project/Krolab/isabella/H3K9me2-Research/ds-analysis"
dir.create(OUT_DIR, showWarnings = FALSE)

PEAK_SCORE_THRESH <- -log10(0.05)
DEG_P_THRESH      <- 0.05

cat("1. Loading data...\n")
peak_gene <- read.table(PEAK_GENE_FILE, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, quote = "", comment.char = "")
cat(sprintf("   peak_gene table: %d rows\n", nrow(peak_gene)))

limma <- read.csv(LIMMA_FILE, stringsAsFactors = FALSE, row.names = 1)
limma$ensembl_id <- rownames(limma)
cat(sprintf("   LIMMA table:     %d genes\n", nrow(limma)))

cat(sprintf("\n2. Defining DMGs: peak_score > %.3f\n", PEAK_SCORE_THRESH))
dmg_set <- unique(peak_gene$gene_id[!is.na(peak_gene$peak_score) &
                                     peak_gene$peak_score > PEAK_SCORE_THRESH])
dmg_set <- dmg_set[!is.na(dmg_set) & dmg_set != ""]
cat(sprintf("   DMG: %d unique genes\n", length(dmg_set)))

cat(sprintf("\n3. Defining DEGs: P.Value < %.2f (full LIMMA)\n", DEG_P_THRESH))
deg_set <- unique(limma$ensembl_id[!is.na(limma$P.Value) &
                                    limma$P.Value < DEG_P_THRESH])
deg_set <- deg_set[!is.na(deg_set) & deg_set != ""]
cat(sprintf("   DEG: %d unique genes\n", length(deg_set)))

overlap    <- intersect(dmg_set, deg_set)
dmg_only   <- setdiff(dmg_set, deg_set)
deg_only   <- setdiff(deg_set, dmg_set)

n_dmg      <- length(dmg_set)
n_deg      <- length(deg_set)
n_overlap  <- length(overlap)
n_dmg_only <- length(dmg_only)
n_deg_only <- length(deg_only)

cat(sprintf("\n4. Overlap counts: DMG-only=%d  overlap=%d  DEG-only=%d\n",
            n_dmg_only, n_overlap, n_deg_only))

coord_in_overlap <- character(0); coord_genes <- character(0)
if (file.exists(COORD_FILE)) {
  coord <- read.csv(COORD_FILE, stringsAsFactors = FALSE)
  coord_genes <- unique(coord$Ensembl_ID)
  coord_in_overlap <- intersect(overlap, coord_genes)
  cat(sprintf("   Coordinated genes: %d, of which %d in overlap\n",
              length(coord_genes), length(coord_in_overlap)))
}

cat("\n5. Attaching gene symbols...\n")
id_to_name <- NULL
if (file.exists(ENSEMBL_MAP)) {
  id_to_name <- read.table(ENSEMBL_MAP, sep = "\t", header = FALSE,
                           stringsAsFactors = FALSE, quote = "", comment.char = "")
  names(id_to_name)[1:2] <- c("ensembl_id", "gene_name")
}
attach_names <- function(ids) {
  df <- data.frame(ensembl_id = ids, stringsAsFactors = FALSE)
  if (!is.null(id_to_name)) {
    df <- merge(df, id_to_name[, c("ensembl_id", "gene_name")],
                by = "ensembl_id", all.x = TRUE, sort = FALSE)
  }
  df
}
write.csv(attach_names(overlap),
          file.path(OUT_DIR, "overlap_DMG_and_DEG.csv"), row.names = FALSE)
write.csv(attach_names(dmg_only),
          file.path(OUT_DIR, "DMG_only.csv"), row.names = FALSE)
write.csv(attach_names(deg_only),
          file.path(OUT_DIR, "DEG_only.csv"), row.names = FALSE)

# Venn plot 
cat("\n6. Making Venn...\n")
make_circle <- function(x0, y0, r, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = x0 + r * cos(theta), y = y0 + r * sin(theta))
}

# Radii scaled by sqrt(count) so areas roughly reflect set size
scale <- 1.8 / sqrt(max(n_dmg, n_deg))
r_dmg <- sqrt(n_dmg) * scale
r_deg <- sqrt(n_deg) * scale

# Position circles so they overlap proportionally
overlap_depth <- min(r_dmg, r_deg) * 1.2
center_sep <- r_dmg + r_deg - overlap_depth
x1 <- -center_sep / 2
x2 <-  center_sep / 2

c1 <- make_circle(x1, 0, r_dmg); c1$set <- "DMG"
c2 <- make_circle(x2, 0, r_deg); c2$set <- "DEG"
circles <- rbind(c1, c2)

# Region-count labels: outside-left, center-overlap, outside-right
overlap_x <- x1 + r_dmg - (overlap_depth / 2)

label_df <- data.frame(
  x     = c(x1 - r_dmg * 0.3,        # DMG-only: centered in left lune
            overlap_x,                # overlap: middle of the lens
            x2 + r_deg * 0.5),        # DEG-only: right portion of DEG circle
  y     = c(0, 0, 0),
  label = c(n_dmg_only, n_overlap, n_deg_only)
)

# Set-name labels: placed well above each circle, offset outward to prevent overlap
label_y <- max(r_dmg, r_deg) + 0.55
set_df <- data.frame(
  x     = c(x1 - r_dmg * 0.6, x2 + r_deg * 0.8),
  y     = c(label_y, label_y),
  label = c(sprintf("DMG (n = %d)", n_dmg),
            sprintf("DEG (n = %d)", n_deg))
)

xlim_pad <- 1.2
ylim_pad <- 1.0
xlim_range <- c(x1 - r_dmg - xlim_pad, x2 + r_deg + xlim_pad)
ylim_range <- c(-max(r_dmg, r_deg) - 0.5, label_y + 0.5)

p <- ggplot() +
  geom_polygon(data = circles,
               aes(x = x, y = y, group = set, fill = set),
               alpha = 0.55, color = "grey25", linewidth = 0.5) +
  geom_text(data = label_df, aes(x = x, y = y, label = label),
            size = 7, fontface = "bold", color = "grey15") +
  geom_text(data = set_df, aes(x = x, y = y, label = label),
            size = 5.5, fontface = "bold") +
  scale_fill_manual(values = c("DMG" = "lightsalmon", "DEG" = "darkolivegreen")) +
  coord_fixed(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(title    = "DMGs vs DEGs",
       subtitle = "Peak raw p < 0.05  |  Expression raw p < 0.05") +
  theme_void() +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 16,
                                     margin = margin(b = 4)),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30",
                                     margin = margin(b = 10)),
        plot.margin   = margin(15, 15, 15, 15))

ggsave(file.path(OUT_DIR, "venn_DMG_vs_DEG.png"),
       p, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "venn_DMG_vs_DEG.pdf"),
       p, width = 8, height = 6, bg = "white")
cat(sprintf("   Saved %s/venn_DMG_vs_DEG.{png,pdf}\n", OUT_DIR))

summary_txt <- paste0(
  "DMG vs DEG Venn Summary\n",
  strrep("-", 50), "\n",
  "DMG: peak_score > 1.301 in peak_gene_DEGs_unfiltered.tsv\n",
  "DEG: P.Value < 0.05 from full limma_deg_results.csv\n",
  "Matched on Ensembl IDs.\n",
  strrep("-", 50), "\n",
  sprintf("DMG:              %d\nDEG:              %d\n", n_dmg, n_deg),
  sprintf("DMG only:         %d\nDEG only:         %d\nOverlap:          %d\n",
          n_dmg_only, n_deg_only, n_overlap)
)
if (length(coord_genes) > 0) {
  summary_txt <- paste0(summary_txt,
    sprintf("Coordinated:      %d  (in overlap: %d)\n",
            length(coord_genes), length(coord_in_overlap)))
}
