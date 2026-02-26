# Load libraries
library(tidyverse)
library(rstatix)
library(ggpubr)
library(patchwork)

# Read data
data <- read.csv("alleleFreq_3_8_11_24_data.csv", row.names = 1)

# Preview
head(data)
dim(data)

# Define groups
high_vir <- data[, 1:12]   # High virulence group (12 samples)
low_vir  <- data[, 13:18]  # Low virulence group (6 samples)

# Results container
results <- data.frame(
  SNP = rownames(data),
  p_value = NA,
  p_value_fdr = NA,
  mean_high = rowMeans(high_vir, na.rm = TRUE),
  mean_low  = rowMeans(low_vir, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Mann‑Whitney U test (two‑sided) for each SNP
for (i in 1:nrow(results)) {
  snp <- results$SNP[i]
  g1 <- as.numeric(high_vir[snp, ]); g1 <- g1[!is.na(g1)]
  g2 <- as.numeric(low_vir[snp, ]);  g2 <- g2[!is.na(g2)]
  
  if (length(g1) > 0 && length(g2) > 0) {
    wt <- wilcox.test(g1, g2, alternative = "two.sided", exact = FALSE)
    results$p_value[i] <- wt$p.value
  }
}

# FDR correction
results$p_value_fdr <- p.adjust(results$p_value, method = "fdr")
results$significant <- results$p_value_fdr < 0.05

# Apply filter: mean_low >= 0.5 & mean_high <= 0.5
filtered <- results %>%
  filter(mean_low >= 0.5, mean_high <= 0.5)

# Significant SNPs after filtering
sig_snps <- filtered %>% filter(significant) %>% arrange(p_value_fdr)

cat("Significant SNPs:\n")
print(sig_snps)

# Save full and significant results
write.csv(results, "all_snp_analysis_results.csv", row.names = FALSE)
write.csv(sig_snps, "significant_snps.csv", row.names = FALSE)

# Prepare data for plot
results <- results %>%
  mutate(
    freq_diff = mean_high - mean_low,
    log10_fdr = -log10(p_value_fdr),
    color_cat = case_when(
      freq_diff > -0.5 ~ "Effect size > -0.5",
      significant     ~ "Significant (FDR < 0.05)",
      TRUE             ~ "Not significant"
    )
  )

# Scatter plot: effect size vs -log10(FDR)
scatter <- ggplot(results, aes(x = freq_diff, y = log10_fdr, color = color_cat)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c(
      "Effect size > -0.5"     = "grey40",
      "Significant (FDR < 0.05)" = "red",
      "Not significant"          = "gray"
    ),
    breaks = c("Significant (FDR < 0.05)", "Effect size > -0.5", "Not significant")
  ) +
  labs(x = "Effect size (VAF_High - VAF_Low)",
       y = "-log10(FDR)", color = "Category") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = max(results$freq_diff, na.rm = TRUE) * 0.8,
           y = -log10(0.05) + 0.1, label = "FDR = 0.05", color = "red", size = 3) +
  annotate("text", x = -0.5, y = max(results$log10_fdr, na.rm = TRUE) * 0.8,
           label = "Effect size = -0.5", color = "red", size = 3, hjust = -0.1)

print(scatter)

# SNPs with freq_diff < -0.5
snps_below <- results %>% filter(freq_diff < -0.5) %>% arrange(desc(freq_diff))
cat("\nSNPs with frequency difference < -0.5:\n")
print(snps_below[, c("SNP", "mean_high", "mean_low", "freq_diff", "p_value_fdr", "significant")])
write.csv(snps_below, "freq_diff_lt_-0.5_snps.csv", row.names = FALSE)

# Boxplot of VAF for these SNPs
wide <- snps_below[, c("mean_high", "mean_low")]
long <- pivot_longer(wide, cols = everything(), names_to = "Virulence", values_to = "VAF")
boxp <- ggplot(long, aes(x = Virulence, y = VAF)) +
  geom_boxplot(width = 0.5, fill = "grey", outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.3, color = "grey20") +
  labs(x = "Virulence", y = "VAF") +
  theme_minimal()

# Combine plots
combined <- scatter + boxp + plot_layout(widths = c(2, 1))
ggsave("dot-box-high-vs-low.png", combined, width = 8, height = 4, dpi = 300)
ggsave("dot-box-high-vs-low.pdf", combined, width = 8, height = 4)

# Summary statistics
cat("\n=== Summary statistics ===\n")
cat("Total SNPs:", nrow(results), "\n")
cat("Significant (FDR < 0.05):", sum(results$significant, na.rm = TRUE), "\n")
cat("freq_diff > -0.5:", sum(results$freq_diff > -0.5, na.rm = TRUE), "\n")
cat("freq_diff <= -0.5:", sum(results$freq_diff <= -0.5, na.rm = TRUE), "\n")

cat("\nBy effect size range:\n")
cat("  > 0      :", sum(results$freq_diff > 0, na.rm = TRUE), "\n")
cat("  0 to -0.5:", sum(results$freq_diff <= 0 & results$freq_diff > -0.5, na.rm = TRUE), "\n")
cat("  <= -0.5  :", sum(results$freq_diff <= -0.5, na.rm = TRUE), "\n")

cat("\nSNPs with effect size > -0.5:\n")
cat("  Total:", nrow(snps_below), "\n")
cat("  Significant:", sum(snps_below$significant, na.rm = TRUE), "\n")
cat("  Non‑significant:", sum(!snps_below$significant, na.rm = TRUE), "\n")

# Top 10 largest freq_diff
cat("\nTop 10 largest freq_diff:\n")
print(results %>% arrange(desc(freq_diff)) %>% head(10) %>%
        select(SNP, mean_high, mean_low, freq_diff, p_value_fdr, significant))

# Bottom 10 smallest freq_diff
cat("\nBottom 10 smallest freq_diff:\n")
print(results %>% arrange(freq_diff) %>% head(10) %>%
        select(SNP, mean_high, mean_low, freq_diff, p_value_fdr, significant))

# Filtered subset: mean_low >= 0.5 & mean_high <= 0.5
cat("\nFiltered SNPs (mean_low >= 0.5 & mean_high <= 0.5):\n")
cat("  Total:", nrow(filtered), "\n")
cat("  Significant:", sum(filtered$significant, na.rm = TRUE), "\n")

# Write summary to file
sink("analysis_summary.txt")
cat("SNP analysis summary\n")
cat("Generated:", Sys.time(), "\n\n")
cat("Total SNPs:", nrow(results), "\n")
cat("Significant (FDR < 0.05):", sum(results$significant, na.rm = TRUE), "\n")
cat("freq_diff > -0.5:", sum(results$freq_diff > -0.5, na.rm = TRUE),
    "(significant:", sum(results$freq_diff > -0.5 & results$significant, na.rm = TRUE), ")\n")
cat("freq_diff <= -0.5:", sum(results$freq_diff <= -0.5, na.rm = TRUE),
    "(significant:", sum(results$freq_diff <= -0.5 & results$significant, na.rm = TRUE), ")\n")
sink()