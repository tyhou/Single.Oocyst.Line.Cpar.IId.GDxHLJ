# Load libraries
library(UpSetR)
library(ggplot2)
library(openxlsx)
library(dplyr)

# Read gene lists
list1 <- scan("list1-3-8-11-24-genes.txt", what = character(), quiet = TRUE)
list2 <- scan("list2-xpehh-genes.txt", what = character(), quiet = TRUE)
list3 <- scan("list3-fst-genes.txt", what = character(), quiet = TRUE)
list4 <- scan("list4-highlyPolymorphic-genes.txt", what = character(), quiet = TRUE)

# Create binary matrix
all_genes <- unique(c(list1, list2, list3, list4))
binary_matrix <- data.frame(
  Gene = all_genes,
  List1 = as.integer(all_genes %in% list1),
  List2 = as.integer(all_genes %in% list2),
  List3 = as.integer(all_genes %in% list3),
  List4 = as.integer(all_genes %in% list4)
)

# Prepare data for UpSetR (genes as row names)
upset_data <- binary_matrix[, -1]
rownames(upset_data) <- binary_matrix$Gene

# Generate UpSet plot
p1 <- upset(upset_data,
            sets = colnames(upset_data),
            nsets = 4,
            order.by = "freq",
            decreasing = TRUE,
            mainbar.y.label = "Number of genes",
            sets.x.label = "Number of genes"
)

# Save plots
png("upset.png", width = 6, height = 4, units = "in", res = 300)
print(p1)
dev.off()

pdf("upset.pdf", width = 6, height = 4)
print(p1)
dev.off()

# Create Excel table with Yes/No indicators
excel_output <- binary_matrix %>%
  mutate(across(starts_with("List"), ~ ifelse(. == 1, "Yes", "")))
write.xlsx(excel_output, "gene_lists_intersection.xlsx")