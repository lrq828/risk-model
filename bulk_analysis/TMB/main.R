###### 获取体细胞突变数据，输出为maf格式#######
# BiocManager::install("PoisonAlien/TCGAmutations")
TCGAmutations::tcga_available()
laml <- TCGAmutations::tcga_load(study = "STAD")
# laml <- TCGAmutations::tcga_load(study = "STAD", source = "Firehose")

###### 计算TMB######
library(maftools)
plotmafSummary(maf = laml,addStat = 'median')

oncoplot(maf = laml, top = 20)

tmb_table_wt_log = tmb(maf = laml)

tmb_table_no_log = tmb(maf = laml,logScale = F)
write.csv(tmb_table_no_log,"tmb_table_no_log.csv")
