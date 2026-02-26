#BiocManager::install("bahlolab/moimix", build_vignettes = TRUE)
library(SeqArray)
setwd("D:/Genome-seq/IId-cross-Single-Oocyst-24-10-07/1_QC_mapping_het_fws/fws/")
#To convert a tabixed gzipped VCF file to the GDS format use the SeqArray package.
seqVCF2GDS("./SOLine_hap_Marker_1065.vcf",
           "./SOLine_hap_Marker_1065.gds")

#To read the GDS file into R, use the seqOpen function and to check that its output matches the orginal VCF file use seqSummary.
my_vcf <-seqOpen("SOLine_hap_Marker_1065.gds")
seqSummary(my_vcf)

# save sample identifiers
sample.id <- seqGetData(my_vcf, "sample.id")
library(moimix)
# get genomic coordinates of all variants
coords <- getCoordinates(my_vcf)
head(coords)

#Estimating the BAF matrix
#The first step in our clonality analysis is to estiamte the B-allele frequencies for each isolate directly from the sequencing depth.
isolate_baf <- bafMatrix(my_vcf)
class(isolate_baf)
str(isolate_baf)
#亲本都是单一的SNV，不混合
plot(isolate_baf, "11730")
plot(isolate_baf, "22971")
#1号小鼠
plot(isolate_baf, "W372")
plot(isolate_baf, "W455")
#24号小鼠
plot(isolate_baf, "W441")
plot(isolate_baf, "W381")
plot(isolate_baf, "W442")
plot(isolate_baf, "W444")
plot(isolate_baf, "W445")
plot(isolate_baf, "W446")
#30号小鼠
plot(isolate_baf, "W388")
plot(isolate_baf, "W464")

#Estimating MOI with Fws
#An alternative way of assessing clonality is to estimte the Fws statistic. An Fws<0.95 is indicative of a clonal infection.
fws_all <- getFws(my_vcf)
hist(fws_all)
# see if our sample that we estimated is multiclonal according to fws
fws_all["W455"] < 0.95
aa <- as.data.frame(fws_all[])
write.table(aa,file = "Fws_SOLine_hap_Marker_1065.txt")
