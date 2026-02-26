library("vcfR")
setwd("D:/Genome-seq/IId-cross-Single-Oocyst-24-10-07/3_AlleleFreq/F1SOSS/")
N <- 14

vcf <- read.vcfR("11SO3SS_hap_Marker_1065.vcf")
chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
ref <- getREF(vcf)
alt <- getALT(vcf)

ad <- extract.gt(vcf, "AD")
ref_split <- masplit(ad, record = 1, sort = 0)
alt_split <- masplit(ad, record = 2, sort = 0)
gt <- extract.gt(vcf, "GT")

datalist <- vector("list", N)

df <- data.frame(CHROM = chrom,
                 POS = pos,
                 REF = ref,
                 ALT = alt,
                 AD_REF.22971 = ref_split[,"22971"],
                 AD_ALT.22971 = alt_split[,"22971"],
                 AD_REF.11730 = ref_split[,"11730"],
                 AD_ALT.11730 = alt_split[,"11730"])

for (i in 1:N) {
  datalist[[i]] <- data.frame(CHROM = chrom,
                              POS = pos,
                              REF = ref,
                              ALT = alt,
                              AD_REF.22971 = ref_split[,"22971"],
                              AD_ALT.22971 = alt_split[,"22971"],
                              AD_REF.11730 = ref_split[,"11730"],
                              AD_ALT.11730 = alt_split[,"11730"],
                              AD_REF.test = ref_split[,i+2],
                              AD_ALT.test = alt_split[,i+2]
  )
}

mask <- which(gt[,"22971"] != "0" &  gt[,"11730"] == "0")
df <- df[mask,]
#构建每个样品的所有画图位点
ref_split_final <- ref_split[mask,]
alt_split_final <- alt_split[mask,]
alt_allele_freq <- alt_split_final/(ref_split_final+alt_split_final)

alt_allele_freq <- as.data.frame(alt_allele_freq)
alt_allele_freq <- cbind(rownames(alt_allele_freq), alt_allele_freq)
library(tidyr)

alt_allele_freq <- separate(data = alt_allele_freq, col = 1, into = c("CHROM", "POS"), sep = "__")
#保存根据reads深度的ALT-allele-freq
library(writexl)
write_xlsx(alt_allele_freq, "alt_allele_freq.xlsx",format_headers = T)
write.table (alt_allele_freq, file ="alt_allele_freq.txt", sep ="\t", row.names =F, col.names =T, quote =F)
alt_allele_freq$POS <- as.numeric(alt_allele_freq$POS)
df <- alt_allele_freq
#把>0.8和<0.2的赋值为1和0，中间值赋值为0.5,效果不好
#df <- as.matrix(alt_allele_freq[,-c(1:2)])
#str(df)
#df[df <= 0.2] <- 0
#df[df >= 0.8] <- 1
#df[df < 0.8 & df > 0.2] <- 0.5
#df <- cbind(alt_allele_freq[,c(1:2)],df)

#计算重组交叉，把allele-freq大于0.6或小于-0.6的认为出现重组交叉
#整理表格
crossover_df1 <- alt_allele_freq[-1065,-1]
crossover_df2 <- alt_allele_freq[-1,-1]
crossover_df1$POS <- as.numeric(crossover_df1$POS)
crossover_df2$POS <- as.numeric(crossover_df2$POS)

crossover <- crossover_df2-crossover_df1
crossover <- cbind(rownames(crossover),rownames(crossover), crossover)
colnames(crossover)[3] <- "length"
colnames(crossover)[2] <- "ID"
crossover <- separate(data = crossover, col = 1, into = c("CHROM", "POS"), sep = "__")
write_xlsx(crossover, "crossover.xlsx",format_headers = T)

library("ggplot2")
library("ggchicklet")
chr <- read.table("Length_info_11730-TGS.txt",header = 1)

picNames <- read.table('new_pic_create.txt') #读取需要新建的图片名称
ID <- read.table('ID.txt') #读取需要的ID

pplist <- vector("list", N)
for (i in 1:N) {
  pplist[[i]] <- ggplot() + 
  geom_point(data = df, aes_string(x= "CHROM", y="POS/1000000", color = ID[i,]), shape = 95, size = 8, alpha = 1) + 
  scale_colour_gradient2(low = "#B71C1C", mid = "#FFCA28", high = "#33691E", midpoint = 0.5) +
  geom_chicklet(data = chr,aes(x= CHROM, y=Length/1000000),width = 0.5,alpha=0,color = 'black') +
  labs(x = '', y = 'Chromosomal Position(Mb)')  + 
  theme_bw() + 
  scale_x_discrete(labels=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8"))+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank (),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) + 
  theme(legend.position = 'none')+
  guides(color = guide_legend(title = FALSE)) 
  
  ggsave(paste("pic_freq_Marker1065/", picNames[i,] ,"_Freq.png", sep =""),plot = pplist[[i]], dpi=300, width=4, height=6, device="png")
  ggsave(paste("pic_freq_Marker1065/", picNames[i,] ,"_Freq.pdf", sep =""),plot = pplist[[i]], dpi=300, width=4, height=6, device="pdf")
  
}

col_GD <- "#33691E"
col_HLJ <- "#B71C1C"

p1 <- ggplot() + 
  geom_chicklet(data = chr,aes(x= CHROM, y=Length/1000000),width = 0.5,alpha=0,color = 'black') +
  geom_point(data = df, aes(x= CHROM, y=POS/1000000), color = col_HLJ, shape = 95, size = 8, alpha = 1) + 
  labs(x = '', y = 'Chromosomal Position(Mb)')  + 
  theme_bw() + 
  scale_x_discrete(labels=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8"))+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank (),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(color = guide_legend(title = FALSE)) 
p1
ggsave("pic_freq_Marker1065/HRB.png",plot = p1, dpi=300, width=4, height=6, device="png")

p2 <- ggplot() + 
  geom_chicklet(data = chr,aes(x= CHROM, y=Length/1000000),width = 0.5,alpha=0,color = 'black') +
  geom_point(data = df, aes(x= CHROM, y=POS/1000000), color = col_GD, shape = 95, size = 8, alpha = 1) + 
  labs(x = '', y = 'Chromosomal Position(Mb)')  + 
  theme_bw() + 
  scale_x_discrete(labels=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8"))+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank (),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) + 
  guides(color = guide_legend(title = FALSE)) 
p2
ggsave("pic_freq_Marker1065/GD.png",plot = p2, dpi=300, width=4, height=6, device="png")
