# Read chromosome length file
chrom_sizes <- read.table("genome.chrom.sizes", header = FALSE, 
                          col.names = c("chrom", "length"))

# Read crossover position data
your_data <- read.table("crossover-15lines-pos-only.txt", header = TRUE)

# Check data structure
cat("Chromosome length file format:\n")
head(chrom_sizes)

cat("\nCrossover data format:\n")
head(your_data)

cat("\nChromosome name comparison:\n")
cat("Chromosomes in chrom.sizes:", unique(chrom_sizes$chrom), "\n")
cat("Chromosomes in crossover data:", unique(your_data$CHROM), "\n")

# Direct matching (both are numbers, just different column names)
your_data$chrom_length <- chrom_sizes$length[match(your_data$CHROM, chrom_sizes$chrom)]

# Check matching result
cat("\nMatching result:\n")
cat("Number of successfully matched sites:", sum(!is.na(your_data$chrom_length)), "/", nrow(your_data), "\n")

if (sum(!is.na(your_data$chrom_length)) == 0) {
  # If still no match, check data types
  cat("\nData type check:\n")
  cat("chrom.sizes$chrom type:", class(chrom_sizes$chrom), "\n")
  cat("your_data$CHROM type:", class(your_data$CHROM), "\n")
  
  # Force conversion to same type
  your_data$CHROM <- as.integer(your_data$CHROM)
  chrom_sizes$chrom <- as.integer(chrom_sizes$chrom)
  
  your_data$chrom_length <- chrom_sizes$length[match(your_data$CHROM, chrom_sizes$chrom)]
  cat("Number of matched sites after type conversion:", sum(!is.na(your_data$chrom_length)), "/", nrow(your_data), "\n")
}

# Remove unmatched rows
your_data <- your_data[!is.na(your_data$chrom_length), ]
cat("Final number of sites for analysis:", nrow(your_data), "\n")



# Function to calculate normalized telomere distance
calculate_telomere_distance <- function(position, chrom_length) {
  # Distance to the nearest telomere
  d_tel <- pmin(position, chrom_length - position)
  # Maximum possible distance (midpoint to telomere)
  max_d <- chrom_length / 2
  # Normalized distance: 0 = telomere, 1 = center
  d_norm <- d_tel / max_d
  return(d_norm)
}

# Calculate telomere distance for each crossover site
your_data$tel_dist <- mapply(calculate_telomere_distance, 
                             your_data$POS, your_data$chrom_length)

# Basic statistics
cat("=== Analysis summary ===\n")
cat("Total number of crossover sites:", nrow(your_data), "\n")
cat("Number of chromosomes involved:", length(unique(your_data$CHROM)), "\n")
cat("Mean normalized telomere distance:", round(mean(your_data$tel_dist), 4), "\n")
cat("Median normalized telomere distance:", round(median(your_data$tel_dist), 4), "\n")
cat("Standard deviation:", round(sd(your_data$tel_dist), 4), "\n")

# Statistical test
t_test <- t.test(your_data$tel_dist, mu = 0.5, alternative = "less")
cat("\n=== Statistical test results ===\n")
cat("t-test - testing whether enriched at ends:\n")
cat("t-value:", round(t_test$statistic, 4), "\n")
cat("p-value:", format.pval(t_test$p.value, digits = 4), "\n")
cat("Observed mean:", round(mean(your_data$tel_dist), 4), "\n")
cat("Expected value under random distribution: 0.5\n")

# Visualization
library(ggplot2)

# Distribution histogram
p1 <- ggplot(your_data, aes(x = tel_dist)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "#C5CAE9", 
                 color = "#5C6BC0", alpha = 0.7) +
  geom_density(color = "#B71C1C", linewidth = 1) +
  geom_vline(xintercept = 0.5, color = "grey20", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = mean(your_data$tel_dist), color = "grey20", linewidth = 1) +
  labs(title = "",
       x = "Normalized distance to telomeres",
       y = "Density") +
  scale_x_continuous(expand = c(0, 0)) +  # x-axis tight
  scale_y_continuous(expand = c(0, 0)) +  # y-axis tight
  theme_classic() 

print(p1)

# Distribution by chromosome
p2 <- ggplot(your_data, aes(x = factor(CHROM), y = tel_dist)) +
  geom_boxplot(fill = "#C5CAE9",color = "#5C6BC0") +
  geom_hline(yintercept = 0.5, color = "grey20", size = 1, linetype = "dashed") +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
  labs(title = "",
       x = "Chromosome",
       y = "Normalized distance to telomeres") +
  theme_classic() 

print(p2)

# Save plots
ggsave("crossover_histogram.png", p1, width = 4, height = 3, dpi = 300)
ggsave("crossover_histogram.pdf", p1, width = 4, height = 3, device = "pdf")
ggsave("crossover_barplot.png", p2, width = 4, height = 3, dpi = 300)
ggsave("crossover_barplot.pdf", p2, width = 4, height = 3, device = "pdf")

# Interpretation of results
cat("\n=== Interpretation of results ===\n")
if (t_test$p.value < 0.05) {
  if (mean(your_data$tel_dist) < 0.5) {
    cat("Statistically significant: crossover sites are significantly enriched at chromosome ends (p < 0.05)\n")
    enrichment_ratio <- (0.5 - mean(your_data$tel_dist)) / 0.5 * 100
    cat("Degree of enrichment toward telomeres:", round(enrichment_ratio, 1), "%\n")
  } else {
    cat("Statistically significant: crossover sites are significantly enriched at chromosome center\n")
  }
} else {
  cat("Not significant: no significant distribution preference found for crossover sites (p > 0.05)\n")
}

# Save results
write.csv(your_data, "crossover_analysis_results.csv", row.names = FALSE)
cat("\nAnalysis complete! Results saved to crossover_analysis_results.csv\n")