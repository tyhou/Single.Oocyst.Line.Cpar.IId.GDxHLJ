plink --vcf SOLine_dip_Marker_1065.vcf --double-id --allow-extra-chr --distance square ibs --out SOLine_distance

# 从 .mibs 生成正确的距离矩阵

# .mibs 是 IBS 共享比例方阵（值越大表示越相似），我们需要将其转换为遗传距离矩阵，公式为：距离 = 1 - 共享比例
awk '{for(i=1;i<=NF;i++) printf "%.6f%s", 1-$i, (i==NF? "\n" : "\t")}' \
    SOLine_distance.mibs > SOLine_distance_matrix.txt

# 创建新环境并安装 quicktree和 r-ape 但后面没用上quicktree
#conda create -n ape_env -c bioconda r-ape -y

# 激活环境
#conda activate ape_env 

# 建树
R --slave --vanilla -e "
  library(ape)
  # 读取 IBS 相似度矩阵
  mibs <- as.matrix(read.table('SOLine_distance.mibs'))
  # 读取样本名：只取第一列（IID）
  sample_names <- read.table('SOLine_distance.mibs.id', header = FALSE)[, 1]
  # 赋予行名和列名
  rownames(mibs) <- sample_names
  colnames(mibs) <- sample_names
  # 转换为遗传距离并建树
  tree <- nj(as.dist(1 - mibs))
  # 输出 Newick 文件
  write.tree(tree, 'SOLine_nj_tree_final.nwk')
"