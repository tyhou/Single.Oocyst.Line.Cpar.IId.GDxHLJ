# Sliding window crossover count plot - CSV version
library(ggplot2)
library(dplyr)
library(scales)   # For axis formatting

# Read CSV data
data <- read.csv("crossover_window_analysis.csv", header = TRUE)

# Check data structure
head(data)
str(data)

# Check column names and rename if needed
# Adjust column names if they differ
colnames(data) <- c("CHROM", "Start", "End", "Crossover_Count")

# Calculate window center position (convert to kb)
data$Window_Center_kb <- (data$Start + data$End) / 2 / 1000
data$Window_Size <- data$End - data$Start

# Classify counts for two-stage plotting
data$count_type <- ifelse(data$Crossover_Count == 0, "zero", "non-zero")

ggplot(data, aes(x = Window_Center_kb, y = Crossover_Count)) +
  # Plot zero-count points first (gray)
  geom_point(data = subset(data, count_type == "zero"), 
             color = "gray", alpha = 0.5, size = 2) +
  # Then plot non-zero points with color gradient
  geom_point(data = subset(data, count_type == "non-zero"), 
             aes(color = Crossover_Count, alpha = Crossover_Count), size = 2.5) +
  scale_color_gradient(low = "blue", high = "red", 
                       name = "Crossover Count\n(non-zero)") +
  scale_alpha_continuous(range = c(0.5, 0.9), guide = "none") +
  scale_x_continuous(labels = comma) +   # Disable scientific notation, use commas
  facet_wrap(~ CHROM, ncol = 1) +
  labs(x = "Genomic Position (kb)", y = "Crossover Count",
       title = "Crossover Count Distribution",
       subtitle = "Gray: Count = 0; Blue-Red: Count > 0") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.line = element_line(color = "black")   # Add axis lines
  )

# Save plots
ggsave("crossover_points.png", width = 10, height = 8, dpi = 300)
ggsave("crossover_points.pdf", width = 10, height = 8, device = "pdf")
