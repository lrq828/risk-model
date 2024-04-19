setwd("E:/GC/bluk/肿瘤微环境/ESTIMATE")
library(utils)
# install.packages("estimate", 
#                  repos="http://r-forge.r-project.org", 
#                  dependencies=TRUE)
library(estimate)
exp <- read.table(file="../genes_exp.txt",header=TRUE, sep="\t",check.names=F,quote="")
head(exp[,1:6],n = 10)
#    GeneSymbol    KD2.1     KD2.2    KD2.3    KD2.4      NC1
## 1      MT-CO1 9531.049 10189.325 8656.926 8873.921 7066.662
## 2      EEF1A1 8108.423  7966.862 7224.281 7783.502 7813.425
## 3       H2-K1 4106.400  4244.788 4031.117 3993.283 5161.887
## 4       H2-D1 4212.534  4433.726 4348.630 4237.329 4455.402
input_file_dir <- '../CIBERSORT/genes_exp.txt'#设置输入的表达数据路径



output_file_dir <- './samples.gct'
output_estimate <- './estimate_score.gct'



filterCommonGenes(input.f = input_file_dir,
                  output.f = output_file_dir,
                  id = "GeneSymbol")
## [1] "Merged dataset includes 9881 genes (531 mismatched)."
estimateScore(input.ds = output_file_dir,
              output.ds= output_estimate, 
              platform="illumina")#注意平台，如果为TCGA或者测序数据则选择illumina
## [1] "1 gene set: StromalSignature  overlap= 136"
## [1] "2 gene set: ImmuneSignature  overlap= 140"




scores <- read.table(output_estimate,skip = 2,header = T)
rownames(scores)=scores[,1]
head(scores[,1:5])

##                        NAME   Description TCGA.D7.5577.01A TCGA.D7.6818.01A TCGA.BR.4280.01A
## StromalScore   StromalScore  StromalScore         122.0569         1165.043        -400.1521
## ImmuneScore     ImmuneScore   ImmuneScore        2063.2019         1518.864        1214.2061
## ESTIMATEScore ESTIMATEScore ESTIMATEScore        2185.2588         2683.907         814.0540
#将表格转置
scores <- data.frame(t(scores[,3:ncol(scores)]))
head(scores)
##                  StromalScore ImmuneScore ESTIMATEScore
## TCGA.D7.5577.01A     122.0569   2063.2019     2185.2588
## TCGA.D7.6818.01A    1165.0429   1518.8643     2683.9072
## TCGA.BR.4280.01A    -400.1521   1214.2061      814.0540
## TCGA.D7.8572.01A     916.3916   1225.5619     2141.9536
## TCGA.VQ.A91Z.01A   -1681.5747   -742.5855    -2424.1602
##TCGA.HU.A4HD.01A    -267.7901    144.7989     -122.9912
scores$sample <- rownames(scores)
scores <- data.frame(scores,row.names = NULL)
write.table(scores,file = 'estimateScore.txt',
            quote = FALSE,sep = "\t",row.names = FALSE)




library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(lemon)
#设置分组信息

scores <- data.table::fread('estimateScore_Group.txt')
#将数据转为长表格，方便绘制分面图
scores2 <- pivot_longer(data = scores,cols = colnames(scores)[1:3],names_to = "typescore",values_to = 'score')
pdf("ESTIMATE_boxplot.pdf",width = 9,height = 7)
ggplot(scores2, aes(x = group, y = score, fill = group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_bw(base_size = 16) + ###去除背景颜色
  theme(panel.grid=element_blank()) +
  facet_rep_wrap(. ~ typescore,scales = 'free',repeat.tick.labels = 'left',ncol = 4)+
  scale_fill_nejm() +
  stat_compare_means(comparisons = combn(unique(scores$group), 2, simplify =FALSE),
                     method = 'wilcox.test')#还可以选择‘wilcox’检验方法
dev.off()






