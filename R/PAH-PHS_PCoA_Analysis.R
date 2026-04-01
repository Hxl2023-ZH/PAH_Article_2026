################################################################################
#                          KEGG PCoA Analysis
################################################################################
# Xinliang Hu
# 202603
# R version 4.2.1

setwd("/path/")

library(dplyr)
library(vegan)
library(ggplot2)
library(pairwiseAdonis)


### caculate relative abudance
counts <- read.table("metadata_merge.txt", header=T, check.names = FALSE, row.names = 1)
metadata <- read.delim('metadata.txt', stringsAsFactors = FALSE,check.names=F)

#加载数据，丰度表格行为样本名，列为feature
otu <- counts[,-1]

#计算Brar_curtis
otu.distance <- vegdist(otu)

#pcoa分析
pcoa <- cmdscale (otu.distance,eig=TRUE) 
pc12 <- pcoa$points[,1:2] 
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2) 

#数据格式转换及数据整合
pc12 <- as.data.frame(pc12)
pc12$Sample <- row.names(pc12)
head(pc12)
p <- ggplot(pc12,aes(x=V1, y=V2))+
  geom_point(size=3)+theme_bw()
p
#colnames(group) <- c("samples","group","Country")
colnames(metadata)

df <- merge(pc12,metadata,by="Sample")
head(df)

col_tmp <- c("#98D7CC", "#FFFFBA", "#C4C0DD", "#FB8C7F", "#8CB8D7")

color=c("#98D7CC", "#FB8C7F", "#C4C0DD", "#8CB8D7")
Shapes <- c(0, 15, 1, 16, 2, 17, 9, 18)

p1<-ggplot(data=df,aes(x=V1,y=V2,
                       color=City,shape=Type))+
  theme_bw()+
  geom_point(size=6)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  #geom_text(aes(label=samples, y=V2+0.03,x=V1+0.03,  vjust=0),size=3.5)+
  #guides(color=guide_legend(title=NULL))+
  labs(x=paste0("PCoA1 ",pc[1],"%"),
       y=paste0("PCoA2 ",pc[2],"%"))+
  scale_color_manual(values = color) +
  scale_shape_manual(values = Shapes) +
  scale_fill_manual(values = c('#00688B', '#CD5B45'))+
  theme(axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20,angle=90),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        legend.position = c(0.8, 0.5),
        panel.grid=element_blank())
  #stat_ellipse(data=df,geom = "polygon",level=0.95,linetype = 2,size=0.5,aes(fill=Group),alpha=0.2,show.legend = T)
p1

#按分组进行PERMANOVA分析
################################ 统计检验 ######################################
# PERMANOVA 默认假设：各组的样本离散度（组内距离）相似
# 离散度
class(otu.distance)
bd <- betadisper(as.dist(otu.distance), df$City) # 方差齐性
anova(bd) 

permutest(bd, permutations = 999)

### 按分组进行PERMANOVA分析
out.adoins0 = adonis2(otu.distance~df$City, distance = 'bray', permutations = 999)

out.adoins1 = adonis2(otu ~ Group, data = metadata, permutations = 999, distance = 'bray')
out.adoins2 = adonis2(otu ~ Type, data = metadata, permutations = 999, distance = 'bray')

out.adoins0
out.adoins1
out.adoins2

### 
otu.dist <- as.dist(as.matrix(otu.distance))
pw1 <- pairwise.adonis2(otu.dist ~ Type, data=df, permutations=999)
pw1

# 依据统计情况添加
p1 + annotate('text', label='PERMANOVA:\n***', x= 0.6, y=0.4, size=5)

ggsave(filename = "Phs+Pah_bray_c_PCoA.png", width = 7, height = 7, units = "in", dpi = 300)
ggsave(filename = "Phs+Pah_bray_c_PCoA.pdf", width = 7, height = 7)






