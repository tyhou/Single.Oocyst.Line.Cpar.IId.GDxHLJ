library(tidyverse)
library(rehh)
library(ggrepel)
library(patchwork)


for(i in 1:8) {
  hap_high_file = paste("phasing/phasing_high_chr", i, ".vcf.gz", sep = "")
  hap_low_file = paste("phasing/phasing_low_chr", i, ".vcf.gz", sep = "")
  hap_all_file = paste("phasing/phasing_all_chr", i, ".vcf.gz", sep = "")
  
  hh_high <- data2haplohh(hap_file = hap_high_file,
                          polarize_vcf = FALSE,
                          vcf_reader = "data.table")
  hh_low <- data2haplohh(hap_file = hap_low_file,
                         polarize_vcf = FALSE,
                         vcf_reader = "data.table")
  hh_all <- data2haplohh(hap_file = hap_all_file,
                         polarize_vcf = FALSE,
                         vcf_reader = "data.table")
  
  scan_high <- scan_hh(hh_high, discard_integration_at_border = FALSE)
  scan_low <- scan_hh(hh_low, discard_integration_at_border = FALSE)
  scan_all <- scan_hh(hh_all, discard_integration_at_border = FALSE)
  
  if (i == 1) {
    wgscan_high <- scan_high
    wgscan_low <- scan_low
    wgscan_all <- scan_all
  } else {
    wgscan_high <- rbind(wgscan_high, scan_high)
    wgscan_low <- rbind(wgscan_low, scan_low)
    wgscan_all <- rbind(wgscan_all, scan_all)
  }
}

wgscan.high.ihs <- ihh2ihs(wgscan_high)
wgscan.low.ihs <- ihh2ihs(wgscan_low)
wgscan.all.ihs <- ihh2ihs(wgscan_all)

write.table(wgscan.high.ihs[["ihs"]], "high.chr.ish.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wgscan.low.ihs[["ihs"]], "low.chr.ish.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wgscan.all.ihs[["ihs"]], "all.chr.ish.txt", sep = "\t", quote = FALSE, row.names = FALSE)

res.xpehh <- ies2xpehh(wgscan_high, wgscan_low, "high", "low")
res.xpehh.2 <- ies2xpehh(wgscan_high, wgscan_all, "high", "all")

write.table(res.xpehh, "high.low.chr.xpehh.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res.xpehh.2, "high.all.chr.xpehh.txt", sep = "\t", quote = FALSE, row.names = FALSE)

if("LOGPVALUE" %in% colnames(res.xpehh)) {
  res.xpehh$LOGPVALUE <- res.xpehh$LOGPVALUE
} else if("PVALUE" %in% colnames(res.xpehh)) {
  res.xpehh$LOGPVALUE <- -log10(res.xpehh$PVALUE)
} else {
  pval_cols <- grep("P|p", colnames(res.xpehh), value = TRUE)
  if(length(pval_cols) > 0) {
    res.xpehh$LOGPVALUE <- -log10(res.xpehh[[pval_cols[1]]])
  } else {
    stop("No p-value column found")
  }
}
res.xpehh$SNP_ID <- paste0("Chr", res.xpehh$CHR, "__", res.xpehh$POSITION)

identify_candidate_regions <- function(xpehh_data, 
                                       p_threshold = 0.05,
                                       merge_distance_kb = 200) {
  logp_threshold <- -log10(p_threshold)
  
  xpehh_data <- xpehh_data %>%
    mutate(
      IS_SIGNIFICANT = LOGPVALUE > logp_threshold,
      SELECTION_DIRECTION = case_when(
        XPEHH_high_low > 0 & IS_SIGNIFICANT ~ "positive",
        XPEHH_high_low < 0 & IS_SIGNIFICANT ~ "negative",
        TRUE ~ "non-significant"
      )
    )
  
  significant_snps <- xpehh_data %>% filter(IS_SIGNIFICANT) %>% arrange(CHR, POSITION)
  
  if (nrow(significant_snps) == 0) {
    return(list(
      xpehh_data = xpehh_data,
      candidate_regions = data.frame(),
      logp_threshold = logp_threshold,
      significant_snps = significant_snps
    ))
  }
  
  merge_distance <- merge_distance_kb * 1000
  candidate_regions <- data.frame()
  current_region <- NULL
  
  for (i in 1:nrow(significant_snps)) {
    snp <- significant_snps[i, ]
    
    if (is.null(current_region)) {
      current_region <- list(
        CHR = snp$CHR,
        START = snp$POSITION,
        END = snp$POSITION,
        N_SNPS = 1,
        MIN_LOGPVALUE = snp$LOGPVALUE,
        MAX_LOGPVALUE = snp$LOGPVALUE,
        MEAN_XPEHH = snp$XPEHH_high_low,
        LEAD_SNP = snp$SNP_ID,
        LEAD_SNP_POS = snp$POSITION,
        LEAD_SNP_LOGPVALUE = snp$LOGPVALUE,
        DIRECTION = ifelse(snp$XPEHH_high_low > 0, "positive", "negative"),
        SNP_LIST = snp$SNP_ID
      )
    } else {
      same_chr <- current_region$CHR == snp$CHR
      within_distance <- snp$POSITION - current_region$END <= merge_distance
      
      if (same_chr && within_distance) {
        current_region$END <- snp$POSITION
        current_region$N_SNPS <- current_region$N_SNPS + 1
        current_region$MIN_LOGPVALUE <- min(current_region$MIN_LOGPVALUE, snp$LOGPVALUE)
        if (snp$LOGPVALUE > current_region$LEAD_SNP_LOGPVALUE) {
          current_region$LEAD_SNP <- snp$SNP_ID
          current_region$LEAD_SNP_POS <- snp$POSITION
          current_region$LEAD_SNP_LOGPVALUE <- snp$LOGPVALUE
        }
        current_region$MAX_LOGPVALUE <- max(current_region$MAX_LOGPVALUE, snp$LOGPVALUE)
        current_region$MEAN_XPEHH <- mean(c(current_region$MEAN_XPEHH, snp$XPEHH_high_low))
        current_region$SNP_LIST <- paste(current_region$SNP_LIST, snp$SNP_ID, sep = ";")
        current_direction <- ifelse(snp$XPEHH_high_low > 0, "positive", "negative")
        if (current_region$DIRECTION != current_direction) {
          current_region$DIRECTION <- "mixed"
        }
      } else {
        candidate_regions <- rbind(candidate_regions, as.data.frame(current_region, stringsAsFactors = FALSE))
        current_region <- list(
          CHR = snp$CHR,
          START = snp$POSITION,
          END = snp$POSITION,
          N_SNPS = 1,
          MIN_LOGPVALUE = snp$LOGPVALUE,
          MAX_LOGPVALUE = snp$LOGPVALUE,
          MEAN_XPEHH = snp$XPEHH_high_low,
          LEAD_SNP = snp$SNP_ID,
          LEAD_SNP_POS = snp$POSITION,
          LEAD_SNP_LOGPVALUE = snp$LOGPVALUE,
          DIRECTION = ifelse(snp$XPEHH_high_low > 0, "positive", "negative"),
          SNP_LIST = snp$SNP_ID
        )
      }
    }
  }
  
  if (!is.null(current_region)) {
    candidate_regions <- rbind(candidate_regions, as.data.frame(current_region, stringsAsFactors = FALSE))
  }
  
  if(nrow(candidate_regions) > 0) {
    candidate_regions <- candidate_regions %>%
      mutate(
        REGION_LENGTH = END - START + 1,
        REGION_LENGTH_KB = round(REGION_LENGTH / 1000, 1),
        REGION_ID = paste0("Region_", 1:n(), "_", CHR),
        REGION_LABEL = paste0(CHR, ":", format(START, scientific = FALSE), "-", 
                              format(END, scientific = FALSE))
      ) %>%
      arrange(CHR, START) %>%
      select(REGION_ID, REGION_LABEL, CHR, START, END, REGION_LENGTH, REGION_LENGTH_KB,
             N_SNPS, DIRECTION, MIN_LOGPVALUE, MAX_LOGPVALUE, MEAN_XPEHH,
             LEAD_SNP, LEAD_SNP_POS, LEAD_SNP_LOGPVALUE, SNP_LIST)
  }
  
  return(list(
    xpehh_data = xpehh_data,
    candidate_regions = candidate_regions,
    logp_threshold = logp_threshold,
    significant_snps = significant_snps
  ))
}

results <- identify_candidate_regions(res.xpehh, p_threshold = 0.05, merge_distance_kb = 200)

write.table(results$xpehh_data, "xpehh_with_significance.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

if (nrow(results$candidate_regions) > 0) {
  write.table(results$candidate_regions, "candidate_selection_regions.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

plot_data <- results$xpehh_data
plot_data$POSITION_PLOT <- plot_data$POSITION
plot_data$POSITION_PLOT <- ifelse(plot_data$CHR %in% "2", plot_data$POSITION_PLOT + 880728, plot_data$POSITION_PLOT)
plot_data$POSITION_PLOT <- ifelse(plot_data$CHR %in% "3", plot_data$POSITION_PLOT + 1873522, plot_data$POSITION_PLOT)
plot_data$POSITION_PLOT <- ifelse(plot_data$CHR %in% "4", plot_data$POSITION_PLOT + 2979129, plot_data$POSITION_PLOT)
plot_data$POSITION_PLOT <- ifelse(plot_data$CHR %in% "5", plot_data$POSITION_PLOT + 4085749, plot_data$POSITION_PLOT)
plot_data$POSITION_PLOT <- ifelse(plot_data$CHR %in% "6", plot_data$POSITION_PLOT + 5180782, plot_data$POSITION_PLOT)
plot_data$POSITION_PLOT <- ifelse(plot_data$CHR %in% "7", plot_data$POSITION_PLOT + 6487972, plot_data$POSITION_PLOT)
plot_data$POSITION_PLOT <- ifelse(plot_data$CHR %in% "8", plot_data$POSITION_PLOT + 7810985, plot_data$POSITION_PLOT)

X_axis <- plot_data %>% group_by(CHR) %>% summarise(center = (max(POSITION_PLOT) + min(POSITION_PLOT)) / 2)

plot_data <- plot_data %>%
  mutate(SIGNED_LOGPVALUE = ifelse(XPEHH_high_low > 0, LOGPVALUE, -LOGPVALUE))

manhattan <- ggplot(plot_data, aes(x = POSITION_PLOT, y = SIGNED_LOGPVALUE)) +
  geom_point(aes(color = as.factor(CHR)), alpha = 0.6, size = 1.0) +
  geom_point(data = subset(plot_data, IS_SIGNIFICANT), 
             aes(fill = SELECTION_DIRECTION), 
             color = "black", size = 2.0, shape = 21, stroke = 0.5) +
  scale_color_manual(values = rep(c("black", "gray40"), 8), guide = "none") +  # hide chromosome legend
  scale_fill_manual(
    name = "Selection type",
    breaks = c("positive", "negative"),
    labels = c("Positive", "Purifying"),
    values = c("positive" = "red", "negative" = "blue", "non-significant" = "white")
  ) +
  scale_x_continuous(
    name = "",
    breaks = X_axis$center,
    labels = paste("Chr", X_axis$CHR),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    name = expression(paste("Signed -log"[10], italic(" P"), " value")),
    breaks = seq(-4, 12, by = 2),
    limits = c(-4, 13),
    expand = expansion(mult = 0.05)
  ) +
  geom_hline(yintercept = c(-results$logp_threshold, results$logp_threshold), 
             color = "black", linetype = "dashed", alpha = 0.7, linewidth = 0.5) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 3, shape = 21))
  )

ggsave("xpehh_manhattan.pdf", plot = manhattan, width = 10, height = 3, dpi = 300)
ggsave("xpehh_manhattan.png", plot = manhattan, width = 10, height = 3, dpi = 300)
