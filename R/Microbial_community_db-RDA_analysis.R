################################################################################
#                 Microbial community db-RDA analysis
################################################################################
# Xinliang Hu
# 202603
# R verison 4.5.2

########################### db-RDA analysis ####################################
setwd("/path/")

library(vegan)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Hmisc)
library(reshape2)

pah <- read.delim("pah_tale.txt", row.names=1, check.names = FALSE)


LMW <- c("NAP","Ace","Acy","Flu","Phe","Ant")
HMW <- c("Flua","Pyr","BaA","Chry","BbF","BkF","BaP","IndP","DahA","BghiP")

pah$LMW_PAH <- rowSums(pah[,LMW])
pah$HMW_PAH <- rowSums(pah[,HMW])
pah$Total_PAH <- pah$LMW_PAH + pah$HMW_PAH

#数据标准化
pah_env <- pah[,c("LMW_PAH","HMW_PAH", "Total_PAH")]
pah_env <- log10(pah_env + 1)

#合并环境变量
env <- read.delim("soil_env_table.txt", row.names=1, check.names = FALSE)
env_all <- cbind(env, pah_env)


#群落数据处理
comm0 <- read.delim("otu_table.txt", row.names=1, check.names = FALSE)
comm <- as.data.frame(t(comm0))
comm_hel <- decostand(comm, method="hellinger")

#构建初始模型
full_model <- capscale(
  comm_hel ~ Total_PAH + LMW_PAH + HMW_PAH + clay + silt + sand + pH + TN + TC + 
    TP + TK + OM + Olsen_P + Olsen_K + NH4 + NO3 + CEC,
  data = env_all,
  distance = "bray"
)

### 自动变量选择
null_model <- capscale(
  comm_hel ~ 1,
  data = env_all,
  distance = "bray"
)
#双向逐步选择（默认）
step_model <- ordistep(
  null_model,
  scope = formula(full_model),
  direction = "forward",
  perm.max = 999
)

step_model
# Call: capscale(formula = comm_hel ~ TP + TK + LMW_PAH + NH4 + HMW_PAH + NO3 + TN,
#               data = env_all, distance = "bray")

# 共线性检测
vif.cca(step_model)
# TP        TK   LMW_PAH       NH4   HMW_PAH       NO3        TN 
# 3.559740  4.011614  5.268660  2.141927  6.016626  3.627092 12.880015 

# 最终db-RDA模型
dbrda_model <- capscale(
  comm_hel ~ LMW_PAH + HMW_PAH + NO3 + TK + TP + NH4, #依据上述结果进行调整
  data = env_all,
  distance = "bray"
)

# 显著性检验 Permutation test
anova(dbrda_model, permutations = 999) # 整体 0.001
anova(dbrda_model, by="terms", permutations=999) #各变量
anova(dbrda_model, by="axis", permutations=999) # 各轴

# 提取解释度
eig_vals <- eigenvals(dbrda_model)
var_explained <- eig_vals / sum(eig_vals)

# 提取绘图数据
site_scores <- scores(dbrda_model, display="sites")
env_scores  <- scores(dbrda_model, display="bp")

site_df <- as.data.frame(site_scores)
env_df  <- as.data.frame(env_scores)

site_df$Sample <- rownames(site_df)
env_df$Variable <- rownames(env_df)

# 添加分组
group <- read.delim("metadata.txt", check.names = FALSE)
site_df <- merge(site_df, group, by="Sample")

# 绘图

p <- ggplot(site_df, aes(CAP1, CAP2)) +
  
  geom_point(aes(color=Group, shape=City),
             size=4) +
  scale_color_manual(values = c("#00688B", "#CD5B45"))+
  geom_segment(
    data=env_df,
    aes(x=0,y=0,xend=CAP1,yend=CAP2),
    arrow=arrow(length=unit(0.3,"cm")),
    colour="black"
  ) +
  
  geom_text(
    data=env_df,
    aes(x=CAP1,y=CAP2,label=Variable),
    size=3
  ) +
  
  theme_bw() +
  
  xlab(paste0("db-RDA1 (", round(var_explained[1]*100,2),"%)")) +
  
  ylab(paste0("db-RDA2 (", round(var_explained[2]*100,2),"%)")) +
  
  theme(
    text=element_text(size=14),
    legend.position="right",
    panel.grid = element_blank()
  )

p

summary(dbrda_model)$cont$importance
RsquareAdj(dbrda_model)

ggsave("microbial_community_RDA_analysis.png", p, width = 5, height = 5, dpi = 300)
ggsave( "microbial_community_RDA_analysis.pdf", p, width = 5, height = 5)

############# 方差分解分析（Variance Partitioning Analysis）#################

pah_vpa <- pah[,c("LMW_PAH","HMW_PAH")]
otu_hel <- comm_hel
soil <- env |> select(NO3, NH4, TK, TP)

pah_vpa <- scale(pah_vpa)
soil <- scale(soil)

#运行VPA
vpa_res <- varpart(otu_hel, pah, soil)
print(vpa_res)

#显著性检验
rda_pah <- rda(otu_hel ~ LMW_PAH + HMW_PAH, data = env_all)
anova(rda_pah, permutations = 999)

rda_soil <- rda(otu_hel ~ NO3 + NH4 + TK + TP, data = env_all)
anova(rda_soil, permutations = 999)


#绘图
vals <- vpa_res$part$indfract$Adj.R.square

vpa_df <- data.frame(
  fraction = c("PAHs", "Shared", "Soil"),
  value = vals[c(1,2,3)]
)
vpa_df

ggplot(vpa_df, aes(x=fraction, y=value)) +
  geom_bar(stat="identity") +
  theme_classic() +
  ylab("Explained variance (Adj R²)") +
  xlab("") +
  theme(
    text = element_text(size=14)
  )


