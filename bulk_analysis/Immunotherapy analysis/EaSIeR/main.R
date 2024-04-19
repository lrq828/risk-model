patient_ICBresponse <- read.csv("TIDE-output(331samples).csv",header = TRUE)
ICBresponse <- patient_ICBresponse
##统计免疫响应患者数
# table(ICBresponse$Responder)
# #  NR   R 
# #  225 106 
# library(tidyverse)
# ICBresponse$Responder <- ifelse(str_detect(ICBresponse$Responder,"False"),"NR","R")
ICBresponse1 <- dplyr::pull(ICBresponse, 3) ##将data.frame中的一列转换为向量
names(ICBresponse1) <- ICBresponse$Patient
ICBresponse1[1:10]

# save(ICBresponse,file = "ICBresponse.Rdata")


################# EaSIeR #######
library("easier")
# load("easier_derived_scores&cell-fraction.Rdata")
groupdat <- read.table(file="group.txt",header=TRUE, sep="\t",check.names=F,quote="")

dat_plot <- data.frame(id = ICBresponse$Patient,
                       t = ICBresponse$TIDE)

##免疫反应预测
easier <- easier_derived_scores
easier$group <- groupdat$clust[match(dat_plot$id,groupdat$samID)] 
library(ggpubr)
p11 <- ggboxplot(easier,x = "group",y = "easier_score",color = "group",add = "jitter") + 
  stat_compare_means()
p11


# load("ICBresponse.Rdata")



##发散条形图
## barplot

# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t > 0,'NR','R'),levels=c('R','NR'))
# 排序
library(dplyr)
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
pA <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  # scale_x_continuous(breaks = seq(0,50,100))+ # X 轴每隔 50 个单位显示一个刻度
  scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('TIDE score') + 
  guides(fill=guide_legend(key.width = 3, key.height = 5, nrow = 2, ncol = 1, byrow = TRUE))+ #显示图例
  theme_prism(border = T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(.95, .95),#图例位置
    legend.justification = c("right", "top")#固定右上角
  )
pA

##堆积条形图

dat_plot$group <-groupdat$clust[match(dat_plot$id,groupdat$samID)] 
group <- c("High-risk","High-risk","Low-risk","Low-risk")
Response <- c("NR","R","NR","R")
# table(dat_plot$group[dat_plot$threshold=="R"])
# table(dat_plot$group[dat_plot$threshold=="NR"])
Percent <- c(0.8571,0.1429,0.5198,0.4802)
num <- c(132,22,92,85)
data <- data.frame(group,Response,Percent,num)

##添加统计数值(卡方检验)
NR <- c(132,92)
R <- c(22,85)
dat <- data.frame(R,NR)
rownames(dat) <- c("High-risk","Low-risk")
chisq.test(dat)
# X-squared = 41.317, df = 1, p-value = 1.295e-10

pB = ggplot( data, aes( x = group, weight = Percent, fill = Response))+
  geom_bar( position = "stack")+xlab("X-squared = 41.317, df = 1, p-value = 1.295e-10")
pC <- pB+theme(panel.grid = element_blank())+theme_bw()+
  scale_fill_manual( values = c('#40E0D0','#FF1493'))
pC

##拼图
library(patchwork)
pD <- (pA+pC+p11)+plot_annotation(tag_levels = 'A')
pD
ggsave(pD,filename = "immune-response.pdf",width = 10,height = 6)




# write.csv(ICBresponse1,"ICBresponse1.csv")
# ICBresponse0 <- read.csv("ICBresponse1.csv",row.names = 1)
# ICBresponse1 <- ICBresponse0$x
# names(ICBresponse1) <- row.names(ICBresponse0)
# rm(ICBresponse0)
# 使用有免疫治疗反应数据评估EaSIeR的预测 --------------------------------------------------
output_eval_with_resp <- assess_immune_response(predictions_immune_response = predictions,
                                                patient_response = ICBresponse1,
                                                RNA_tpm = RNA_tpm,
                                                TMB_values = TMB,
                                                easier_with_TMB = "weighted_average",
                                                weight_penalty = 0.5)

# 通过下列方式调出图片并保存
# inspect output
p1 <- output_eval_with_resp[[1]]
p2 <- output_eval_with_resp[[2]]
p3 <- output_eval_with_resp[[3]]


# 通过Biomarker解释对免疫疗法的反应
output_biomarkers <- explore_biomarkers(pathways = pathway_activities,
                                        immunecells = cell_fractions,
                                        tfs = tf_activities,
                                        lrpairs = lrpair_weights,
                                        ccpairs = ccpair_scores,
                                        cancer_type = 'STAD',
                                        patient_label = ICBresponse1)

save(output_biomarkers,file = "output_biomarkers.Rdata")

load( "output_biomarkers.Rdata")
p4 <- output_biomarkers[[1]];p4
p5 <- output_biomarkers[[2]];p5
p6 <- output_biomarkers[[3]];p6
p7 <- output_biomarkers[[4]];p7
p8 <- output_biomarkers[[5]];p8
p9 <- output_biomarkers[[6]];p9



# 免疫细胞亚型的数据也拿到了，一起放进来

# rm(list = ls())
library(dplyr)
library(tidyr)
##加载分组数据
# load("cmoicova.Rdata")
# groupdat <- cmoic.ova$clust.res
# groupdat$clust <- ifelse(groupdat$clust==1,"C1","C2")

##免疫细胞类型


head(cell_fractions[,1:4])
cell_f <- cell_fractions %>%as.data.frame() %>% 
  mutate(.,samID=rownames(.)) %>% 
  dplyr::select(.,-Other)

cell_group <- cell_f
cell_group$group <-groupdat$clust[match(dat_plot$id,groupdat$samID)] 

#宽转长——gather
cell_dat <- tidyr::gather(data = cell_group,key = "Cell",value = "Fractions",
                          B:DC)
colnames(cell_dat)
##分组柱形图
library(ggplot2)
library(ggpubr)
#  分组柱形图
# https://blog.csdn.net/zhouhucheng00/article/details/106368179
p0 <- ggplot(cell_dat, aes(Cell, Fractions)) +
  geom_boxplot(aes(color = group))+
  scale_color_manual(values = c("#00AFBB", "red")) +
  stat_compare_means(aes(group = group), label = "p.signif")
p10 <- p0+theme_classic()
ggsave(p10,filename = "immune-type.pdf",width = 15,height = 5)


library(patchwork)
p_all <- (p5+p8)/(p6+p7)
p <- p_all+plot_annotation(tag_levels = 'A');p
ggsave(p,filename = "immune-subtype.pdf",width = 25,height =20)
