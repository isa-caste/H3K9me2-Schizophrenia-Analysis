# load packages 
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(ggpubr)    
  library(grid)
})
 # directory set up
PROJECT_ROOT <- "/N/project/Krolab/isabella/H3K9me2-Research/ds-analysis"
OUT_DIR <- PROJECT_ROOT
CORR_CSV <- file.path(PROJECT_ROOT, "per_sample_correlations.csv")
corr <- read_csv(CORR_CSV, show_col_types = FALSE)

# function to format p-values
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(sprintf("p = %.2e", p))
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
}

############################################
# auto detect column names
###########################################
group_col <- intersect(c("group","Group","sample_group","SampleGroup",
                         "condition","Condition"), names(corr))[1]
r_col <- intersect(c("correlation","pearson_r","r","negative_pearson",
                     "neg_pearson_r","neg_r","pearson"), names(corr))[1]
if (is.na(group_col) || is.na(r_col)) {
  stop("Couldn't autodetect columns in per_sample_correlations.csv.\n",
       "Columns present: ", paste(names(corr), collapse = ", "),
       "\nEdit group_col / r_col manually.")
}
message("  group column: ", group_col, " | value column: ", r_col)

corr$group_clean <- ifelse(
  grepl("cont|npc|ctrl", corr[[group_col]], ignore.case = TRUE),
  "Control", "Case")

############################################
# Stat tests
###########################################
ctrl_vals <- corr[[r_col]][corr$group_clean == "Control"]
case_vals <- corr[[r_col]][corr$group_clean == "Case"]

# Shapiro-Wilk normality tests
sw_ctrl <- shapiro.test(ctrl_vals)
sw_case <- shapiro.test(case_vals)
# Wilcoxon rank-sum test (non-parametric)
wx <- wilcox.test(ctrl_vals, case_vals, exact = FALSE)
# T-test (parametric, for reference)
tt <- t.test(ctrl_vals, case_vals, var.equal = FALSE)

#summary
summary_txt <- paste0(
  "per-sample Pearson r: Control vs Case\n",
  strrep("-", 46), "\n",
  sprintf("n Control = %d, n Case = %d\n", length(ctrl_vals), length(case_vals)),
  sprintf("median r: Control = %.4f, Case = %.4f\n",
          median(ctrl_vals), median(case_vals)),
  sprintf("mean r:   Control = %.4f, Case = %.4f\n",
          mean(ctrl_vals), mean(case_vals)),
  "\n",
  sprintf("Shapiro-Wilk (Control): W = %.4f, %s\n",
          sw_ctrl$statistic, fmt_p(sw_ctrl$p.value)),
  sprintf("Shapiro-Wilk (Case):    W = %.4f, %s\n",
          sw_case$statistic, fmt_p(sw_case$p.value)),
  "  (Shapiro p < 0.05 means the group is NOT normally distributed\n",
  "   -> Wilcoxon is the right test. If both >= 0.05, t-test is fine too.)\n",
  "\n",
  sprintf("Wilcoxon rank-sum: W = %.2f, %s\n",
          wx$statistic, fmt_p(wx$p.value)),
  sprintf("Welch's t-test (reference): t = %.3f, df = %.2f, %s\n",
          tt$statistic, tt$parameter, fmt_p(tt$p.value)),
  "\n",
  sprintf(">> caption: %s (Wilcoxon rank-sum)\n",
          fmt_p(wx$p.value))
)
cat(summary_txt)
writeLines(summary_txt, file.path(OUT_DIR, "wk_stats_summary.txt"))

# create boxplot
corr_flipped <- corr
corr_flipped$correlation_flipped <- -corr[[r_col]]

boxplot_flipped <- ggplot(corr_flipped, aes(x = group_clean, y = correlation_flipped, fill = group_clean)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.7, color = "grey30") +
  scale_fill_manual(values = c(Control = "cornflowerblue",  Case = "tomato3")) +
  labs(x = "Sample Group",
       y = "Negative Pearson Correlation (-r)",
       title = "Pearson Correlation by Group") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(OUT_DIR, "boxplot_per_sample_correlations.png"), boxplot_flipped,
       width = 7, height = 5, dpi = 300)
