setwd("E:/GC/bluk/肿瘤微环境/CIBERSORT")
#install.packages('e1071')
library(ggplot2)
library(reshape2)
library(ggpubr)
library(e1071)
library(dplyr)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library("limma")
#读取输入文件，并对输入文件整理
# setwd("E:/GC/bluk/免疫景观分析/Cibersort")            #设置工作目录

source('Cibersort.R')

# 设置分析依赖的基础表达文件
# 每类免疫细胞的标志性基因及其表达
# 基因名字为Gene symbol
LM22.file <- "LM22.txt"
# 1. Cibersort
TCGA_exp.file <- "genes_exp.txt"

TCGA_TME.results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 100, QN = F)  
# perm置换次数=1000
# QN如果是芯片设置为T，如果是测序就设置为F

write.csv(TCGA_TME.results, "CIBERSORT_Results.csv")


result <- data.table::fread('CIBERSORT_Results.csv')

library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(result,key=CIBERSORT,value = Composition,-c(Group,Sample)) # 提前准备好group信息

pdf("CIBERSORT-boxplot.pdf",width = 10,height = 7)
ggboxplot(b, x = "CIBERSORT", y = "Composition",
          fill = "Group", palette = "lancet")+
  stat_compare_means(aes(group = Group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()
