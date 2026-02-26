# PCA analysis script
# Input: eigenvec and eigenval from PLINK, population info file
# Output: PCA plots with eigenvalue tables

# Load required libraries
library(ggplot2)
library(ggforce)   # for facet_zoom
library(dplyr)

# Read eigenvec and eigenval
eigvec <- read.table("SOLine_dip_Marker_1065.eigenvec", header = FALSE, stringsAsFactors = FALSE)
eigval <- read.table("SOLine_dip_Marker_1065.eigenval", header = FALSE)

# Compute variance explained by each PC
percentage <- eigval$V1 / sum(eigval$V1) * 100
pcs <- paste0("PC", 1:nrow(eigval))

# Create eigenvalue summary table
eigval_df <- data.frame(PCs = pcs,
                        variance = eigval$V1,
                        proportion = percentage,
                        stringsAsFactors = FALSE)
write.table(eigval_df, file = "plink.eigenvalue.xls", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Read population assignment file
poptable <- read.table("SOLine_pca.pop.txt", header = TRUE, comment.char = "")
# Assume the 4th column contains group labels
pop <- unique(poptable[, 4])

# Prepare PCA data (first 5 eigenvectors)
pca.data <- eigvec[, 2:6]
colnames(pca.data) <- c("vcf_id", "PC1", "PC2", "PC3", "PC4")
pca.data <- merge(pca.data, poptable, by = "vcf_id")

# Set factor levels for groups (custom order)
group_levels <- c("HLJ", "GD", "F1", "Line1", "Line2", "Line3", "Line4", "Line6", "Line7",
                  "Line8", "Line10", "Line11", "Line16", "Line19", "Line24",
                  "Line30", "Line31", "Line32")
pca.data$Group <- factor(pca.data$Group, levels = group_levels)

# Define colour palette for groups
c3 <- c("#E05A5A", "#BAD563", "#FFCD33", "#9BC2E6", "#F4B183", "#C15811", "#757575", "#F06292",
        "#BA124A", "#BA8CDC", "#7030A0", "#BF7F3F", "#0964AF", "#774F27", "#063F6E",
        "#68A141", "#02192C", "#385723")

# Extract variance percentages for PC1 and PC2
pc1_perc <- round(percentage[1], 1)
pc2_perc <- round(percentage[2], 1)

# Create main scatter plot of PC1 vs PC2
p12 <- ggplot(pca.data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(alpha = 0.8, size = 4) +
  theme_bw() +
  scale_color_manual(values = c3) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  labs(x = paste0("PC1 (", pc1_perc, "%)"),
       y = paste0("PC2 (", pc2_perc, "%)"))

# Add zoomed view for a specific region
p12_zoom <- p12 + facet_zoom(xlim = c(-0.12, -0.02), ylim = c(-0.05, 0.05),
                             zoom.size = 0.5)

# Save plots
ggsave("SOLine_PCA-PC12-zoom.pdf", plot = p12_zoom, device = "pdf",
       dpi = 300, width = 8, height = 4)
ggsave("SOLine_PCA-PC12-zoom.png", plot = p12_zoom, device = "png",
       dpi = 300, width = 8, height = 4)