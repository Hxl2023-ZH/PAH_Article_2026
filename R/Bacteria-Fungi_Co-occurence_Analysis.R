################################################################################
#                       Bacteria-Fungil Co-occurence Analysis
################################################################################
# Xinliang Hu
# 202603
# R version 4.2.1

library(Hmisc)
library(tidyverse)
library(pheatmap)
library(igraph)
library(RColorBrewer)
library(ggraph)
library(scales)

#------------------------------ 加载数据 ---------------------------------------
# 输入文件路径
setwd('/path/')

bac_file  <- "otu_tax_table.txt"
fun_file  <- "ITS_otu_table.txt"

# 读取
bac0 <- read.table(bac_file, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
bac0[1:5, 1:5]

# filter
keep <- rowSums(bac0 > 0) >= ncol(bac0)*0.2
bac0 <- bac0[keep, ]
dim(bac0)
# [1] 583  24

fun0 <- read.table(fun_file, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
fun0[1:5, 1:5]

# filter
keep <- NULL
keep <- rowSums(fun0 > 0) >= ncol(fun0)*0.2
fun0 <- fun0[keep, ]
dim(fun0)
# [1] 176  24

metadata <- read.table("metadata.txt", 
                       sep = '\t', header = TRUE, check.names = FALSE)

# 提取分组数据
metadata_Poluted <- metadata %>% filter(Group == "clean")
sp <- metadata_Poluted$Sample

bac <- bac0[,sp]
fun <- fun0[,sp]

# 对齐样本
common_samples <- intersect(colnames(bac), colnames(fun))
bac <- bac[, common_samples]
fun <- fun[, common_samples]


#-------------计算细菌与真菌之间的 Spearman 相关（使用 rcorr）------------------

# 转置为 sample × taxa 格式（rcorr 要变量在列）
bac_t <- t(bac)
fun_t <- t(fun)

res_list <- list()

for (b in colnames(bac_t)) {
  for (f in colnames(fun_t)) {
    x <- bac_t[, b]
    y <- fun_t[, f]
    # 跳过零方差
    if (sd(x)==0 | sd(y)==0) {
      rho <- NA; p <- NA
    } else {
      ct <- rcorr(cbind(x, y), type="spearman")
      rho <- ct$r[1,2]
      p   <- ct$P[1,2]
    }
    res_list[[length(res_list)+1]] <-
      data.frame(Bacteria=b, Fungi=f, Rho=rho, Pvalue=p)
  }
}

corr_df <- bind_rows(res_list)
corr_df$FDR <- p.adjust(corr_df$Pvalue, method="BH")

# 输出
write.table(corr_df, "Bacteria_Fungi_correlation_table.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


#---------------------- 筛选显著的共现关系（边）--------------------------------

cut_rho <- 0.7
cut_fdr <- 0.05

edge_df <- corr_df %>%
  filter(!is.na(Rho)) %>%
  filter(abs(Rho) >= cut_rho & FDR <= cut_fdr)


write.table(edge_df, "Bacteria_Fungi_edges_filtered.txt",
            sep="\t", quote=FALSE, row.names=FALSE)




#--------------------------- 构建共现网络 --------------------------------------

# network edge list
edges <- edge_df %>% 
  mutate(weight = Rho,
         type = ifelse(Rho>0, "positive", "negative"))

# create node list
nodes_bac <- data.frame(name = unique(edges$Bacteria), kingdom="Bacteria")
nodes_fun <- data.frame(name = unique(edges$Fungi),   kingdom="Fungi")
nodes <- bind_rows(nodes_bac, nodes_fun)

# build igraph network
g <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)

E(g)$Rho <- as.numeric(E(g)$Rho)

# ------------------ Louvain 模块划分（子网络分析）-----------------------------
E(g)$weight <- abs(E(g)$Rho)
cl <- cluster_louvain(g)
V(g)$module <- cl$membership

write.table(data.frame(Node=V(g)$name, Kingdom=V(g)$kingdom,
                       Module=V(g)$module),
            "Bacteria_Fungi_node_modules.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


#-------------------------- 网络拓扑参数 --------------------------------------
V(g)$degree       <- degree(g)
V(g)$betweenness  <- betweenness(g)
V(g)$closeness    <- closeness(g)
V(g)$eigenvector  <- eigen_centrality(g)$vector

centrality_df <- data.frame(
  Node=V(g)$name,
  Kingdom=V(g)$kingdom,
  Module=V(g)$module,
  Degree=V(g)$degree,
  Betweenness=V(g)$betweenness,
  Closeness=V(g)$closeness,
  Eigenvector=V(g)$eigenvector
)

write.table(centrality_df, "Bacteria_Fungi_network_centrality.txt",
            sep="\t", quote=FALSE, row.names=FALSE)



# -------------------------- 共现网络图（ggraph） ------------------------------
# color schemes

# ---- 定义形状 ----
shape_map <- c(
  "Bacteria" = 21,  # 圆形
  "Fungi" = 24      # 三角形
)

# ---- 定义模块颜色 ----

#module_colors <- brewer.pal(max(V(g)$module), "Set3") # 根据情况选择

module_colors <- colorRampPalette(brewer.pal(12, "Set3"))(max(V(g)$module))

V(g)$module_color <- module_colors[V(g)$module]

# ---- 绘图 ----

#  ayout = "stress"

p <- ggraph(g, layout = "stress") +
  
  # 边（相关性）
  geom_edge_link(
    aes(color = Rho, width = abs(Rho)),
    alpha = 0.7
  ) +
  scale_edge_color_gradient2(
    low = "#00688B", mid = "white", high = "#CD5B45",
    midpoint = 0, name = "Correlation"
  ) +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  
  # 节点（颜色=模块；形状=Kingdom）
  geom_node_point(
    aes(
      fill = module_color,
      shape = kingdom
    ),
    size = 3,
    stroke = 0.8,
    color = "black"
  ) +
  scale_shape_manual(values = shape_map, name = "Kingdom") +
  scale_fill_identity() +  # 模块颜色不需要图例
  
  # 节点标签
  
  theme_graph(base_family = 'Helvetica') +
  ggtitle("Bacteria–Fungi Co-occurrence Network")

p

ggsave("Cooccurrence_network_with_modules_stress.pdf", p, width = 9, height = 7)
#ggsave("Cooccurrence_network_with_modules_graphopt.pdf", p, width = 9, height = 7)

#———————————————————————————— 文件保存 -----------------------------------------
# 保存边文件
edges_out <- data.frame(
  source = edge_df$Bacteria,
  target = edge_df$Fungi,
  weight = edge_df$Rho,
  pvalue = edge_df$Pvalue,
  FDR = edge_df$FDR
)

write.table(edges_out,
            file = "cytoscape_edges.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


# 保存节点文件
node_df <- data.frame(
  name = V(g)$name,
  kingdom = V(g)$kingdom,
  module = V(g)$module,
  degree = degree(g),
  betweenness = betweenness(g, normalized = TRUE),
  eigen = eigen_centrality(g)$vector
)

write.table(node_df,
            file = "cytoscape_nodes.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


#----------------------------- 子网络比较 --------------------------------------

modules <- split(V(g)$name, V(g)$module)

module_stats <- lapply(names(modules), function(mod){
  
  module_nodes <- modules[[mod]]
  g_sub <- induced_subgraph(g, vids = module_nodes)
  
  data.frame(
    module = as.numeric(mod),
    n_nodes = gorder(g_sub),
    n_edges = gsize(g_sub),
    density = edge_density(g_sub),
    mean_abs_cor = mean(abs(E(g_sub)$Rho)),
    mean_degree = mean(V(g_sub)$degree),
    mean_betweenness = mean(V(g_sub)$betweenness),
    mean_eigenvector = mean(V(g_sub)$eigenvector)
  )
})

module_stats_df <- bind_rows(module_stats)

write.table(module_stats_df,
            "Bac_Fungi_polluted_module_comparison.txt",
            sep="\t", quote=FALSE, row.names = FALSE)



















