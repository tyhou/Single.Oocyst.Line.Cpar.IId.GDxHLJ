# Load required libraries
library(tidyverse)
library(writexl)
library(gridExtra)

# Read nucleotide diversity (Pi) data for low and high virulence groups
dt_low <- read.delim("keep.low_1000win_100step.windowed.pi", sep = "\t", header = TRUE, check.names = FALSE)
dt_high <- read.delim("keep.high_1000win_100step.windowed.pi", sep = "\t", header = TRUE, check.names = FALSE)

# Read FST between high and low virulence groups
fst_high_low <- read.delim("fst_1000win_100step_virluence.windowed.weir.fst", sep = "\t", header = TRUE, check.names = FALSE)

# Merge datasets
dt2 <- merge(dt_low, dt_high, by = c("CHROM", "BIN_START", "BIN_END"), all = TRUE)
dt3 <- merge(dt2, fst_high_low, by = c("CHROM", "BIN_START", "BIN_END"), all = TRUE)

# Rename Pi columns for clarity
names(dt3)[names(dt3) == "PI.x"] <- "PI.low"
names(dt3)[names(dt3) == "PI.y"] <- "PI.high"
dt3[is.na(dt3)] <- 0

# Select relevant columns
df3 <- dt3[, c("CHROM", "BIN_START", "PI.low", "PI.high", "N_VARIANTS", "WEIGHTED_FST")]

# Compute ratio of Pi (high / low) and its log2
df3$ROD <- df3$PI.high / df3$PI.low
df3$log2ROD <- log2(df3$ROD)
df3 <- df3[!is.infinite(df3$log2ROD), ]

# Define thresholds (75th percentile for log2ROD, 25th for others)
vline1 <- quantile(df3$log2ROD, 0.75)      # upper quartile
vline2 <- quantile(df3$log2ROD, 0.25)      # lower quartile
hline <- quantile(df3$WEIGHTED_FST, 0.75)  # upper quartile of FST

# Regions potentially under positive selection in high virulence group:
# high FST and reduced diversity (log2ROD < vline2)
select_gene_high <- df3[df3$log2ROD < vline2 & df3$WEIGHTED_FST > hline, ]

# Regions where high virulence group shows increased diversity (for completeness)
select_gene_low <- df3[df3$log2ROD > vline1 & df3$WEIGHTED_FST > hline, ]

# Save merged data and selected regions to Excel
xl_list <- list(
  merge_pi_Fst = df3,
  select_gene_low = select_gene_low,
  select_gene_high = select_gene_high
)
write_xlsx(xl_list, "merge_pi_Fst-1000win-100step.xlsx", format_headers = FALSE)

# Main scatter plot: log2(ROD) vs FST
p_main <- ggplot(df3, aes(x = log2ROD, y = WEIGHTED_FST)) +
  # Highlight the top-left quadrant (candidate positive selection in high virulence)
  annotate("rect", xmin = -Inf, xmax = vline2, ymin = hline, ymax = Inf,
           fill = "red", alpha = 0.2) +
  geom_point(size = 2, alpha = 0.8) +
  # Threshold lines
  geom_hline(yintercept = hline, color = "red", linewidth = 1, linetype = "dotted") +
  geom_vline(xintercept = vline1, color = "red", linewidth = 1, linetype = "dotted") +
  geom_vline(xintercept = vline2, color = "red", linewidth = 1, linetype = "dotted") +
  labs(x = expression(log[2]("Pi"[High] / "Pi"[Low])), y = "Fst") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_text(size = rel(0.8)),
    axis.text = element_text(size = rel(0.8), color = "black"),
    plot.margin = unit(c(0, 0, 0, 0), "inches")
  )

# Top marginal density plot (log2ROD)
p_top <- ggplot(df3, aes(x = log2ROD)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#dedede", colour = "#dedede", binwidth = 0.08) +
  theme_classic() +
  labs(y = "Density") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Right marginal density plot (FST)
p_right <- ggplot(df3, aes(x = WEIGHTED_FST)) +
  geom_histogram(aes(y = after_stat(density)), fill = "#dedede", colour = "#dedede", binwidth = 0.01) +
  theme_classic() +
  labs(y = "Density") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_flip()  # flip to match orientation

# Empty placeholder for top-right corner
empty <- ggplot() + geom_point(aes(1, 1), colour = "white") +
  theme_void()

# Arrange plots: top marginal, empty, main, right marginal
p_combined <- grid.arrange(
  p_top, empty, p_main, p_right,
  ncol = 2, nrow = 2,
  widths = c(4, 1), heights = c(1, 4)
)

# Save output
ggsave("ROD_Fst_low_high_1000win_100step_density.pdf", plot = p_combined,
       device = "pdf", dpi = 300, width = 6, height = 6)
ggsave("ROD_Fst_low_high_1000win_100step_density.png", plot = p_combined,
       device = "png", dpi = 300, width = 6, height = 6)