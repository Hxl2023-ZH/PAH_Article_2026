################################################################################
#                      Bac-Fungi lm Analysis
################################################################################
# Xinliang Hu
# 202603
# R verison 4.2.1

setwd("/path/")

library(ggpmisc)
library(ggplot2)
library(ggpubr)

df = read.table("Bac_ITS_diversity_lm_table.txt", header = TRUE, sep = "\t", 
                row.names = 1, stringsAsFactors = FALSE)
head(df)


#Remove confidence region (se = FALSE)
# Extend the regression lines: fullrange = TRUE
draw1 <- ggplot(df, aes(x = ITS_Pielou_e, y = ITS_Pielou_e))+
  geom_point(aes(color = Group, shape = Group)) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "top",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20))+
  #geom_rug(aes(color =Group)) +
  geom_smooth(aes(color = Group),
              method = lm, 
              se = TRUE, 
              fullrange = TRUE)+
  scale_color_manual(values = c('#00688B', '#CD5B45'))+
  ggpubr::stat_cor(method = "spearman", aes(color = Group), label.x = 1)+
  scale_y_continuous(expand = c(0,0.1))+
  #scale_x_continuous(limits = c(0, 15), expand = c(0.01,0.01))+
  labs(x = "Pielou_e (Bacteria)", y = 'Pielou_e (Fungi)\n')

draw1

ggsave(draw1, filename = "Fungi_bacteriome_richness_lm.png", width = 7, 
       height = 7, units = "in", dpi = 300)

ggsave(draw1, filename = "Fungi_bacteriome_Pielou_lm.pdf", width = 7, height = 7)



