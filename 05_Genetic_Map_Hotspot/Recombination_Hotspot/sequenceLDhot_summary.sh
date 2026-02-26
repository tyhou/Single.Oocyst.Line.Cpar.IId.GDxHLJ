

sbatch sequenceLDhot.sh

cd output_W20_Q75
R
source("../HotspotSummary.R")
HotspotSummary("./chr1_sequenceLDhot.sum",plot = T,output = T,ofile = "./chr1_hotspot.txt",title = "chr1")
HotspotSummary("./chr2_sequenceLDhot.sum",plot = T,output = T,ofile = "./chr2_hotspot.txt",title = "chr2")
HotspotSummary("./chr3_sequenceLDhot.sum",plot = T,output = T,ofile = "./chr3_hotspot.txt",title = "chr3")
HotspotSummary("./chr4_sequenceLDhot.sum",plot = T,output = T,ofile = "./chr4_hotspot.txt",title = "chr4")
HotspotSummary("./chr5_sequenceLDhot.sum",plot = T,output = T,ofile = "./chr5_hotspot.txt",title = "chr5")
HotspotSummary("./chr6_sequenceLDhot.sum",plot = T,output = T,ofile = "./chr6_hotspot.txt",title = "chr6")
HotspotSummary("./chr7_sequenceLDhot.sum",plot = T,output = T,ofile = "./chr7_hotspot.txt",title = "chr7")
HotspotSummary("./chr8_sequenceLDhot.sum",plot = T,output = T,ofile = "./chr8_hotspot.txt",title = "chr8")
q()