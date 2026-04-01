################################################################################
#                       MAGs Corralation Analysis
################################################################################
# Xinliang Hu
# 202603
# R version: 4.2.1

setwd("/path/")

library(Hmisc)
library(tidyverse)
library(pheatmap)
library(igraph)
library(RColorBrewer)
library(ggraph)
library(scales)

##############################    Load File   ##################################
set.seed(2025)

abund <- read.table("MAGs_4_city_relative_abundance_final.txt", 
                    sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

metadata <- read.table("metadata.txt", 
                       sep = '\t', header = TRUE, check.names = FALSE)

metadata_Poluted <- metadata %>% filter(Group == "Polluted")
sp <- metadata_Poluted$SampleID

abund0 <- abund[,sp]

# only keep MAGs that present at least in 4 samples(Polluted group)
abund0_1 <- abund0
abund0_1[abund0_1 > 0] <- 1
abund0_1[1:5,]
abund0_1$Total <- apply(abund0_1, 1, FUN = sum)

abund0$Total <- abund0_1$Total
abund0 <- abund0 %>% filter(Total > 3)
abund0 <- abund0[,-13]
dim(abund0) # Polluted [1] 77 12

all(colnames(abund0) %in% metadata_Poluted$SampleID)

############################# Corralation Analysis #############################

# cor caculator
cor_mat <- rcorr(t(abund0), type = "spearman")

head(cor_mat)


#将矩阵转换为数据框
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

df1 <- flattenCorrMatrix(cor_mat$r,cor_mat$P)


# Benjamini-Hochberg correction
df1$FDR <- p.adjust(df1$p, method = "BH")
head(df1)
#df2 <- df1[df1$cor > 0.3,]

#df2 <- df1 %>% filter(cor > 0.3 | cor < -0.3) %>% filter(p < 0.05) %>% filter(FDR < 0.05)
df2 <- df1 %>% filter(p < 0.05) %>% filter(FDR < 0.05)
head(df2)

write.table(df2,file="Cor_polluted_filter.txt",
            sep = "\t",row.names=FALSE, col.names=T, quote = FALSE)

# 所有涉及的 MAG
mags <- sort(unique(c(df2$row, df2$column)))

# 创建空矩阵
mat_df <- matrix(NA, nrow = length(mags), ncol = length(mags))
rownames(mat_df) <- mags
colnames(mat_df) <- mags

# 填充对称相关系数
for (i in 1:nrow(df2)) {
  r1 <- df2$row[i]
  r2 <- df2$column[i]
  val <- df2$cor[i]
  mat_df[r1, r2] <- val
  mat_df[r2, r1] <- val
}

# 对角线 = 1
diag(mat_df) <- 1

mat_df[is.na(mat_df)] <- 0

pdf("MAGs_heatmap_polluted15.pdf", width = 32, height = 32)
pheatmap(mat_df,
         color = colorRampPalette(c("#00688B", "white", "#CD5B45"))(100),
         clustering_method = "average",
         border_color = "black",
         main = "MAGs Correlation Heatmap",
         na_col = "grey90",
         #filename = "MAG_correlation_heatmap.png",
         width = 12, height = 12)
dev.off()

############################ Network Analysis ##################################

edge_df <- df2 %>%
  filter(abs(cor) >= 0.7, FDR < 0.05) %>%
  select(row, column, cor)

g <- graph_from_data_frame(edge_df, directed = FALSE)
g

E(g)$color <- ifelse(E(g)$cor > 0, "#CD5B45", "#00688B")
E(g)$weight <- abs(E(g)$cor)
E(g)$width  <- rescale(E(g)$weight, to = c(0.3, 2.5))  # 缩放线宽

# plot
p <- ggraph(g, layout = "stress") +    # stress 布局更均匀
  geom_edge_link(aes(color = cor, width = width),
                 alpha = 0.8,
                 show.legend = TRUE) +
  scale_edge_color_gradient2(
    low = "#00688B", mid = "white", high = "#CD5B45",
    midpoint = 0,
    name = "Correlation"
  ) +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  geom_node_point(shape = 21, fill = "gray", size = 5, color = "black") +
  geom_node_text(aes(label = name),
                 repel = TRUE,               # 关键！防止标签重叠
                 size = 3.5,
                 point.padding = unit(0.2, "lines")) +
  theme_graph(base_family = "sans") +
  ggtitle("MAG Correlation Network")

p
ggsave("MAGs_polluted_network_raw.pdf", p, width = 7, height = 7)

# raw plot
plot(g,
     vertex.size = 8,
     vertex.label.cex = 0.6,
     edge.curved = 0.1,
     layout = layout_with_fr(g),
     main = "MAGs Correlation Network")

########################### module analysis ####################################
# 使用 Louvain 聚类
comm <- cluster_louvain(g)

# 每个节点所属子网络（模块）
V(g)$module <- comm$membership

sizes(comm) #查看模块数量与大小：

n_mod <- max(V(g)$module)
module_colors <- brewer.pal(min(n_mod, 12), "Set3")  # 最多 12 个颜色
if (n_mod > 12) module_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_mod)

V(g)$color <- module_colors[V(g)$module]

# 保存Louvain模块信息
modules <- split(V(g)$name, V(g)$module)

write.table(
  data.frame(
    module = rep(names(modules), sapply(modules, length)),
    MAG = unlist(modules)
  ),
  "MGAs_module_polluted_list.txt",
  sep="\t", quote=FALSE, row.names = FALSE
)

# raw plot
plot(g,
     vertex.size = 10,
     vertex.label.cex = 0.7,
     layout = layout_with_fr,
     main = "MAG Correlation Network (Louvain Modules)")


# 边颜色（正红 / 负蓝）
E(g)$color <- ifelse(E(g)$cor > 0, "red", "blue")
E(g)$width <- rescale(abs(E(g)$cor), to = c(0.3, 2.5))

# 绘制美化后的总网络图
p2 <- ggraph(g, layout = "stress") +
  geom_edge_link(aes(color = cor, width = width),
                 alpha = 0.7, show.legend = TRUE) +
  scale_edge_color_gradient2(
    low = "#00688B", mid = "white", high = "#CD5B45",
    midpoint = 0,
    name = "Edge Correlation"
  ) +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  
  geom_node_point(aes(color = factor(module)),
                  size = 6,
                  alpha = 0.9,
                  show.legend = TRUE) +
  scale_color_manual(values = module_colors, name = "Module") +
  
  geom_node_text(aes(label = name),
                 repel = TRUE,
                 size = 3.5,
                 max.overlaps = 20,
                 point.padding = unit(0.15, "lines")) +
  
  theme_graph(base_family = "sans") +
  ggtitle("MAG Correlation Network (Louvain Modules)")

p2

ggsave("MAGs_polluted_network_Modules.pdf", p2, width = 9, height = 7)

###### 提取最大模块（节点数最多
largest_mod <- which.max(sizes(comm))

largest_mod <- 4
sub_nodes <- V(g)$name[V(g)$module == largest_mod]

g_sub <- induced_subgraph(g, vids = sub_nodes)

plot(g_sub,
     vertex.size = 12,
     vertex.label.cex = 0.8,
     edge.width = abs(E(g_sub)$cor) * 5,
     edge.color = E(g_sub)$color,
     layout = layout_with_fr(g_sub),
     main = paste("Subnetwork - Module", largest_mod))

col_tmp <- c("#98D7CC", "#FFFFBA", "#C4C0DD", "#FB8C7F", "#8CB8D7")
### ggraph plot
p_sub <- ggraph(g_sub, layout = "stress") +
  geom_edge_link(aes(color = cor, width = width),
                 alpha = 0.8, show.legend = TRUE) +
  scale_edge_color_gradient2(
    low = "#00688B", mid = "white", high = "#CD5B45",
    midpoint = 0,
    name = "Correlation"
  ) +
  geom_node_point(size = 9, color = "#FB8C7F", fill = "gold") +
  geom_node_text(aes(label = name),
                 repel = TRUE,
                 size = 6,
                 point.padding = unit(0.2, "lines")) +
  theme_graph(base_family = 'Helvetica') +
  ggtitle(paste("Major Subnetwork ( Module", largest_mod, ")"))
p_sub
ggsave("MAGs_polluted_network_Modules_M04.pdf", p_sub, width = 7, height = 7)

#### 保存完整网络的边表/节点表
edge_df_final <- data.frame(
  from = ends(g, E(g))[,1],
  to   = ends(g, E(g))[,2],
  cor  = E(g)$cor,
  module_from = V(g)$module[match(ends(g, E(g))[,1], V(g)$name)],
  module_to   = V(g)$module[match(ends(g, E(g))[,2], V(g)$name)]
)

write.table(edge_df_final,
            "network_edges.txt",
            sep="\t", quote=FALSE, row.names = FALSE)


node_df0 <- data.frame(
  MAG    = V(g)$name,
  module = V(g)$module,
  degree = degree(g)
)

write.table(node_df0,
            "MAG_polluted_network_nodes.txt",
            sep="\t", quote=FALSE, row.names = FALSE)

########################### 计算全局网络中心性 #################################
# Degree
V(g)$degree <- degree(g)

# Closeness
V(g)$closeness <- closeness(g, normalized = TRUE)

# Betweenness
V(g)$betweenness <- betweenness(g, normalized = TRUE)

# Eigenvector Centrality
ev <- eigen_centrality(g)$vector
V(g)$eigenvector <- ev


centrality_df <- data.frame(
  MAG        = V(g)$name,
  module     = V(g)$module,
  degree     = V(g)$degree,
  closeness  = V(g)$closeness,
  betweenness= V(g)$betweenness,
  eigenvector= V(g)$eigenvector
)

write.table(centrality_df,
            "MAG_polluted_node_centrality.txt",
            sep="\t", quote=FALSE, row.names = FALSE)

############################## 网络模块比较分析 ################################

modules <- split(V(g)$name, V(g)$module)

module_stats <- lapply(names(modules), function(mod){
  
  module_nodes <- modules[[mod]]
  g_sub <- induced_subgraph(g, vids = module_nodes)
  
  data.frame(
    module = as.numeric(mod),
    n_nodes = gorder(g_sub),
    n_edges = gsize(g_sub),
    density = edge_density(g_sub),
    mean_abs_cor = mean(abs(E(g_sub)$cor)),
    mean_degree = mean(V(g_sub)$degree),
    mean_betweenness = mean(V(g_sub)$betweenness),
    mean_eigenvector = mean(V(g_sub)$eigenvector)
  )
})

module_stats_df <- bind_rows(module_stats)

write.table(module_stats_df,
            "MAG_polluted_module_comparison.txt",
            sep="\t", quote=FALSE, row.names = FALSE)

# 模块规模（节点/边数）柱状图
p_size <- ggplot(module_stats_df, aes(x=factor(module), y=n_nodes)) +
  geom_bar(stat="identity", fill="#00688B", color = "black") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15)
  ) +
  labs(x="Module", y="Number of Nodes",
       title="Node Count per Module")
p_size

ggsave("MAG_polluted_module_node_count.pdf", p_size, width=7, height=7)

# 模块密度比较
p_density <- ggplot(module_stats_df, aes(x=factor(module), y=density)) +
  geom_bar(stat="identity", fill="#00688B", color = "black") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15)
  ) +
  labs(x="Module", y="Network Density",
       title="Network Density per Module")
p_density

ggsave("MAG_polluted_module_density.pdf",p_density, width=8, height=5)

# 模块对比分析：哪个模块最重要？根据模块规模、密度、中心性加权得到“模块重要性”
module_stats_df$importance_score <- 
  scale(module_stats_df$n_nodes) +
  scale(module_stats_df$density) +
  scale(module_stats_df$mean_eigenvector)

module_stats_df <- module_stats_df %>% arrange(desc(importance_score))

write.table(module_stats_df,
            "MAG_polluted_module_importance_ranking.txt",
            sep="\t", quote=FALSE, row.names = FALSE)

################################## Hub MAG #####################################

cent_df <- data.frame(
  MAG = V(g)$name,
  module = V(g)$module,
  degree = V(g)$degree,
  closeness = V(g)$closeness,
  betweenness = V(g)$betweenness,
  eigenvector = V(g)$eigenvector
)

# Top 10% cutoff
deg_cutoff  <- quantile(cent_df$degree, 0.90)
bet_cutoff  <- quantile(cent_df$betweenness, 0.90)
eig_cutoff  <- quantile(cent_df$eigenvector, 0.90)

hub_df <- cent_df %>%
  filter(
    degree >= deg_cutoff |
      betweenness >= bet_cutoff |
      eigenvector >= eig_cutoff
  )

write.table(hub_df,
            "MAG_polluted_HUB_MAG_list.txt",
            sep="\t", quote=FALSE, row.names = FALSE)

# 模块内 Hub MAG
hub_by_module <- hub_df %>% arrange(module, desc(degree))

write.table(hub_by_module,
            "MAG_network_results/HUB_MAG_by_module.txt",
            sep="\t", quote=FALSE, row.names = FALSE)

# 分析不同模块中中心性是否存在差异
p_mod_box <- ggplot(cent_df, aes(x=factor(module), y=degree)) +
  geom_boxplot(fill="#00688B") +
  geom_point(data=hub_df, aes(x=factor(module), y=degree),
             color="#CD5B45", size=3) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15)
  ) +
  labs(title="Degree Distribution by Module\n(Hub MAG Marked in Red)",
       x="Module", y="Degree")
p_mod_box

ggsave("MAG_polluted_module_degree_boxplot_with_HUB.pdf",
       p_mod_box, width=8, height=5)

col_tmp <- c("#98D7CC", "#FFFFBA", "#C4C0DD", "#FB8C7F", "#8CB8D7")
# Hub MAG 气泡图
hub_long <- hub_df %>%
  pivot_longer(cols=c(degree, closeness, betweenness, eigenvector),
               names_to="centrality",
               values_to="value")

p_bubble <- ggplot(hub_long,
                   aes(x=centrality, y=MAG, size=value, color=factor(module))) +
  geom_point(alpha=0.8) +
  scale_color_manual(values = c("#98D7CC", "#FFFFBA", "#C4C0DD", "#FB8C7F")) +
  scale_size_continuous(name="Centrality") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 15))+
  labs(title="Hub MAG Centrality Overview",
       x="Centrality Type", y="MAG")

p_bubble

ggsave("MAG_polluted_HUB_MAG_centrality_bubble.pdf",
       p_bubble, width=7, height=7)














