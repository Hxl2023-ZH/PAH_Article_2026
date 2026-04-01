################################################################################
#                           KEGG Annotation Analysis
################################################################################
# Xinliang Hu
# 202603
# R version 4.2.1

setwd('/path/')

library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)

library(clusterProfiler)
library(MicrobiomeProfiler)

cts <- as.matrix(read.csv('merged_4city_ko_Counts.txt',sep="\t",row.names="KOID"))
cts[1:5, 1:6]

coldata <- read.delim('metadata.txt', row.names=1)
head(coldata)

coldata$Group <- factor(coldata$Group)

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# Chongqing
coldata_CQ <- coldata %>% filter(City == 'CQ')
CQ_selt <- rownames(coldata_CQ)
cts_CQ0 <- cts[,CQ_selt]
head(cts_CQ0)

# JiNan
coldata_JN <- coldata %>% filter(City == 'JN')
JN_selt <- rownames(coldata_JN)
cts_JN <- cts[,JN_selt]
head(cts_JN)

# HeiFei
coldata_HF <- coldata %>% filter(City == 'HF')
HF_selt <- rownames(coldata_HF)
cts_HF <- cts[,HF_selt]

# HZ
coldata_HZ <- coldata %>% filter(City == 'HZ')
HZ_selt <- rownames(coldata_HZ)
cts_HZ <- cts[,HZ_selt]

####################### KO富集分析

res <- NULL
res1 <- NULL
res1_select <- NULL
res1_up <<- NULL
res1_down <- NULL

KO_enrich_analysis <- function(mycts, city, outdir = ".", mycoldata) {
  
  dds <- DESeqDataSetFromMatrix(countData = mycts,
                                colData = mycoldata,
                                design = ~ Group)
  
  keep <- rowSums(counts(dds)) >= 10 # keep only rows that have at least 10 reads total
  dds <- dds[keep,]
  dds
  
  # Setting the factor levels
  dds$Group <- factor(dds$Group, levels = c("polluted","clean"))
  
  dds1 <- DESeq(dds, 
                fitType = 'mean', 
                minReplicatesForReplace = 7, 
                parallel = FALSE)
  
  #注意，需将 treat 在前，control 在后，意为 treat 相较于 control 中哪些基因上调/下调
  res <<- results(dds1, contrast=c("Group","polluted","clean"))
  res
  
  #输出表格至本地
  res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  outfile_res1 <- file.path(outdir, paste0("KO_gene_", city, ".DESeq2.txt"))
  write.table(res1, outfile_res1, col.names = NA, sep = '\t', quote = FALSE)
  
  ##筛选差异基因
  #首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
  res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  
  #log2FC≥1 & padj<0.01 标识 Enriched，代表显著富集的基因
  #log2FC≤-1 & padj<0.01 标识 Deleted，代表显著降低的基因
  #其余标识 none，代表非差异的基因
  res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'Enriched'
  res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'Deleted'
  res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'
  
  res1 <<- res1
  
  #输出选择的差异基因总表
  res1_select <<- subset(res1, sig %in% c('Enriched', 'Deleted'))
  
  outfile_res1_sel <- file.path(outdir, paste0("KO_", city, ".DESeq2.select.txt"))
  write.table(outfile_res1, file = outfile_res1_sel, sep = '\t', col.names = NA, quote = FALSE)
  
  #根据 up 和 down 分开输出
  res1_up <<- subset(res1, sig == 'Enriched')
  res1_down <<- subset(res1, sig == 'Deleted')
  
  outfile_res1_up <- file.path(outdir, paste0("KO_", city, ".DESeq2.Enriched.txt"))
  write.table(res1_up, file = outfile_res1_up, sep = '\t', col.names = NA, quote = FALSE)
  
  outfile_res1_down <- file.path(outdir, paste0("KO_", city, ".DESeq2.Deleted.txt"))
  write.table(res1_down, file = outfile_res1_down, sep = '\t', col.names = NA, quote = FALSE)
  
}


KO_enrich_analysis(mycts = cts_HZ, city = "HZ", outdir = "HZ/", mycoldata = coldata_HZ)

# add 'KO' Colum
res1$KO <- rownames(res1)

# 按‘log2FoldChange’排序提取注释数据
topKOs <- res1 %>% 
  filter(sig == "Enriched") %>% 
  arrange(desc(abs(log2FoldChange))) %>%
  slice(1:10)

# 按‘padj’进行排序提取注释数据
topKOs <- res1 %>% 
  filter(sig == "Enriched") %>% 
  arrange(padj) %>%
  slice(1:10)

downKOs <- res1 %>% 
  filter(sig == "Deleted") %>% 
  arrange(padj) %>%
  slice(1:10)



### ggplot2 差异火山图

res11 <- res1 %>% filter(!is.na(padj))

#默认情况下，横轴展示 log2FoldChange，纵轴展示 -log10 转化后的 padj
p <- ggplot(data = res11, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 2) + #绘制散点图
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
  xlim(-20, 20) + ylim(0, 200)   #定义刻度边界

p

# 添加标签
p1 <- p +
  geom_text_repel(
    data = topKOs,
    aes(x = log2FoldChange, y = -log10(padj), label = KO),
    size = 3.5,
    box.padding = 0.4,
    show.legend = FALSE
  ) +
  geom_text_repel(
    data = downKOs,
    aes(x = log2FoldChange, y = -log10(padj), label = KO),
    size = 3.5,
    box.padding = 0.4,
    show.legend = FALSE
  ) 
p1

ggsave(p1, filename = 'HZ/KO_HZ_volcano_plot_labels.pdf', width = 7, height = 7)

########################## KEGG Enrich Analysis ################################

###### enriched
head(res1_up)
outdir <- "HZ"
type_KO <- "Enriched"
city <- "HZ"


res1_type <- res1_up

res1_type$KOid <- rownames(res1_type)

ekegg <- enrichKEGG(
  gene         = res1_type$KOid,
  organism     = "ko",     # 宏基因组使用 ko
  keyType      = "kegg",
  use_internal_data = F,# KO 基因 ID
  pvalueCutoff = 0.01
)

head(ekegg)


file1 <- file.path(outdir, paste0("KEGG_Enrichment_of_", type_KO, "_", city, ".csv"))
write.csv(as.data.frame(ekegg), file1, row.names = FALSE)


file2 <- file.path(outdir, paste0("KEGG_Enrichment_of_", type_KO, "_", city, "_dotplot.pdf"))
pdf(file2, width = 7, height = 7)
dotplot(ekegg) + 
  ggtitle("KEGG Enrichment of Enriched KOs")+
  scale_color_gradient(low = "#CD5B45", high = "#00688B")

dev.off()

file3 <- file.path(outdir, paste0("KEGG_Enrichment_of_", type_KO, "_", city, "_cbarplot.pdf"))
pdf(file3, width = 7, height = 5)
barplot(ekegg, showCategory = 20)
dev.off()


library(enrichplot)
file4 <- file.path(outdir, paste0("KEGG_Enrichment_of_", type_KO, "_", city, "_cnetplot.pdf"))
pdf(file4, width = 12, height = 12)
cnetplot(ekegg, foldChange = setNames(res1_select$log2FoldChange, res1_select$KOid))
dev.off()


################################# Pathway view #################################
library(pathview)

top_pathways <- ekegg$ID[1:10]

fc = res1_select$log2FoldChange
names(fc) = rownames(res1_select)

pathview(
  gene.data  = fc,
  pathway.id = "00624", 
  species    = "ko", 
  out.suffix = "Polycyclic_aromatic",
  kegg.native=TRUE
)





