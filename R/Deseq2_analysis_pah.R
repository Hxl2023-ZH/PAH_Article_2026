################################################################################
#                           DESeq2 Analysis
################################################################################
# Xinliang Hu
# 2025-11-10

setwd('/path/')

library(DESeq2)
library(ggplot2)

cts <- as.matrix(read.csv('MAG006_JN_gene_read_counts.txt',sep="\t",row.names="gene_id"))
cts[1:5, 1:6]

coldata <- read.delim('JN_metadata.txt', row.names=1)
head(coldata)

coldata$Group <- factor(coldata$Group)

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)] # reorder?
all(rownames(coldata) == colnames(cts))


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group)
dds

keep <- rowSums(counts(dds)) >= 10 # keep only rows that have at least 10 reads total
dds <- dds[keep,]

# Setting the factor levels
dds$Group <- factor(dds$Group, levels = c("polluted","clean"))

dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

#注意，需将 treat 在前，control 在后，意为 treat 相较于 control 中哪些基因上调/下调
res <- results(dds1, contrast=c("Group","polluted","clean"))

res


#输出表格至本地
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res1, 'MAG006_gene_JN.DESeq2.txt', col.names = NA, sep = '\t', quote = FALSE)


##筛选差异基因
#首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

#log2FC≥1 & padj<0.01 标识 Enriched，代表显著富集的基因
#log2FC≤-1 & padj<0.01 标识 Deleted，代表显著降低的基因
#其余标识 none，代表非差异的基因
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'Enriched'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'Deleted'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'

#输出选择的差异基因总表
res1_select <- subset(res1, sig %in% c('Enriched', 'Deleted'))
write.table(res1_select, file = 'MAG006_gene_JN.DESeq2.select.txt', sep = '\t', col.names = NA, quote = FALSE)

#根据 up 和 down 分开输出
res1_up <- subset(res1, sig == 'Enriched')
res1_down <- subset(res1, sig == 'Deleted')

write.table(res1_up, file = 'MAG006_gene_JN.DESeq2.Enriched.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'MAG006_gene_JN.DESeq2.Deleted.txt', sep = '\t', col.names = NA, quote = FALSE)


##ggplot2 差异火山图

res1$KO <- rownames(res1)

# 按‘padj’进行排序提取注释数据
topKOs <- res1 %>% 
  filter(sig == "Enriched") %>% 
  arrange(padj) %>%
  slice(1:10)

downKOs <- res1 %>% 
  filter(sig == "Deleted") %>% 
  arrange(padj) %>%
  slice(1:10)


#默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
p <- ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 3) +  #绘制散点图
  scale_color_manual(values = c('#CD5B45', 'gray', '#00688B'), 
                     limits = c('Enriched', 'none', 'Deleted')) +  #自定义点的颜色
  labs(x = 'Log2 Fold Change', 
       y = '-log10 (adjust p-value)', 
       title = 'Polluted vs Clean', 
       color = '') +  #坐标轴标题
  theme(plot.title = element_text(hjust = 0.5, size = 20), 
        panel.grid = element_blank(), #背景色、网格线、图例等主题修改
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),
        legend.position = c(0.9,0.9),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  geom_vline(xintercept = c(-1, 1), lty = 3, color = 'black') +  #添加阈值线
  geom_hline(yintercept = 2, lty = 3, color = 'black') +
  xlim(-10, 10) + ylim(0, 20)  #定义刻度边界

p

# 添加标签
p1 <- p +
  geom_text_repel(
    data = topKOs,
    aes(x = log2FoldChange, y = -log10(padj), label = KO),
    size = 3.5,
    box.padding = 0.4,
    show.legend = FALSE,
    max.overlaps = 10
  ) +
  geom_text_repel(
    data = downKOs,
    aes(x = log2FoldChange, y = -log10(padj), label = KO),
    size = 3.5,
    box.padding = 0.4,
    show.legend = FALSE
  ) 
p1

ggsave(p1, filename = 'MAG006_gene_JN_volcano_plot_label.pdf', width = 7, height = 7)




