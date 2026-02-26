# 加载必要的包
library(ggplot2)
library(scales)
setwd("D:/Genome-seq/IId-cross-Single-Oocyst-24-10-07/1_QC_mapping_het_fws/Het-cov/")
# 假设你的数据已经存储在名为df的数据框中
# df 包含四列：ID, Heterozygosity, Coverage, Category

# 创建示例数据（如果你没有实际数据，可以用这个测试）
# df <- data.frame(
#   ID = c(paste0("Sample_O_", 1:11), paste0("Sample_S_", 1:3)),
#   Heterozygosity = c(runif(11, 0.2, 0.8), runif(3, 0.1, 0.3)),
#   Coverage = c(runif(11, 30, 90), runif(3, 85, 98)),
#   Category = rep(c("Oocyst", "Sporozoite"), times = c(11, 3))
# )
df <- read.table(file = "het-cov.txt",header = T, row.names = NULL)
# 绘制组合图
p <- ggplot(df, aes(x = ID)) +
  # 1. 绘制柱状图 (左边Y轴: Heterozygosity)
  geom_col(aes(y = Heterozygosity, fill = Category), width = 0.7, alpha = 0.8) +
  # 2. 绘制折线图和点 (右边Y轴: Coverage)
  geom_line(aes(y = Coverage / 100, group = 1), color = "#E57373", linewidth = 1) + 
  geom_point(aes(y = Coverage / 100), color = "#E57373", size = 2, shape = 19) +
  
  # 3. 设置双坐标轴
  scale_y_continuous(
    name = "Heterozygosity", # 左边Y轴名称
    limits = c(0, 1), # 左边Y轴范围
    sec.axis = sec_axis(~. * 100, 
                        name = "Average coverage (%)", 
                        labels = function(x) paste0(x, "%")) # 添加百分号
  ) +
  
  # 4. 设置柱状图的填充颜色
  scale_fill_manual(values = c("Single_oocyst" = "#FFD456", "Single_sporozoite" = "#7F88CB")) +
  
  # 5. 设置主题和标签
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y.left = element_text(color = "black"),
    axis.title.y.left = element_text(color = "black", face = "bold"),
    axis.text.y.right = element_text(color = "#E57373"),
    axis.title.y.right = element_text(color = "#E57373", face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey40"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "",
    x = "",
    fill = "",
    caption = ""
  )

# 显示图形
print(p)

# 如果需要保存图片
ggsave("heterozygosity_coverage_plot.png", p, device = "png", width = 6, height = 4, dpi = 300)
ggsave("heterozygosity_coverage_plot.pdf", p, device = "pdf", width = 6, height = 4, dpi = 300)
