library(stringr)
library(ggplot2)
library(ggpubr)
library(patchwork)


result <- data.table::fread('TIDE-output.csv')
colnames(result)
#根据行名新增分组列：
result$Risk <- ifelse(
  str_sub(result$Patient,1,1)=='L','Low-risk','High-risk'
)
#转换为因子，调整顺序：
result$Risk <- factor(result$Risk,levels = c('Low-risk','High-risk'))
#小提琴图展示结果：
#1.TIDE小提琴图：
my_comparisons <- list( c("Low-risk", "High-risk")) #添加比较分组
p1 <- ggviolin(result, x = 'Risk', y = 'TIDE', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p1
#method可选：t.test、wilcox.test、anova、kruskal.test



#2.Dysfunction小提琴图：
p2 <- ggviolin(result, x = 'Risk', y = 'Dysfunction', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p2



#3.Exclusion小提琴图：
p3 <- ggviolin(result, x = 'Risk', y = 'Exclusion', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 't.test')
p3


#4.MSI小提琴图：
colnames(result)[6]
colnames(result)[6] <- c('MSI') #简化一下列名
p4 <- ggviolin(result, x = 'Risk', y = 'MSI', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", bracket.size=0.5, tip.length = 0.02, method = 'anova')
p4

p <- p1+p2+p3+p4
pdf("TIDE.pdf")
p1
dev.off()
pdf("ALL.pdf")
p
dev.off()
pdf("noMSI.pdf",width = 30,height = 10)
p1+p2+p3
dev.off()


