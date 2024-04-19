gene <- read.delim('gene_bluk.txt',row.names = 1)
miRNA <- read.delim('miRNA_bluk.txt',row.names = 1)
library(Hmisc) #载入包
rcorr_test_spearman <- rcorr(as.matrix(miRNA),as.matrix(gene),type = 'pearson')#这里的编写规则和前两个不一样，数据需要转化为矩阵as.matrix()
rcorr_test_spearman_R <- rcorr_test_spearman$r[1:1881,1882:2201]
rcorr_test_spearman_P <- rcorr_test_spearman$P[1:1881,1882:2201]#后来琢磨了一下，我只选择结果里面我需要的数据即可
write.csv(rcorr_test_spearman_R,'R.csv')#导出数据
write.csv(rcorr_test_spearman_P,'P.csv')
