library(parallel)
library(pRRophetic)
library(ggplot2)
suppressWarnings(suppressMessages(future::plan("multiprocess",workers = 8)))
rm(list = ls())
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
data(cgp2016ExprRma)
possibleDrugs2016 <- unique(drugData2016$Drug.name)
possibleDrugs2016
#用system.time来返回计算所需时间
head(possibleDrugs2016)
expdata <- read.delim("genes_exp.txt")
rownames(expdata)=expdata[,1]
expdata=expdata[,-1]
expdata <- as.matrix(expdata)

all.results <- data.frame(colnames(expdata))
names(all.results) <- "Risk"
library(pRRophetic)
# data(bortezomibData)
for (i in 218:length(possibleDrugs2016)) {
  result=pRRopheticPredict(
    testMatrix=expdata,
    drug=possibleDrugs2016[i],
    tissueType="all",
    batchCorrect="eb",
    selection=1,
    dataset="cgp2016")
    
  result <- as.data.frame(result)
  names(result) <- possibleDrugs2016[i]
  all.results <- cbind(all.results,result)
}
write.csv(all.results,"allDrugssenstivity.csv")

#######25,41,48,52,54,58,89,107,180,183,197,199,213,217



#############
system.time({
  cl <- makeCluster(3)
  results <- parLapply(cl,possibleDrugs2016[1:3],
                     function(x){
                       library(pRRophetic)
                       data(bortezomibData)
                       predictedPtype=pRRopheticPredict(
                         testMatrix=exprDataBortezomib,
                         drug=x,
                         tissueType="all",
                         batchCorrect="eb",
                         selection=1,
                         dataset="cgp2016")
                       return(predictedPtype)
                     })#Lapply.的并行版本
  res.df <- do.call('rbind',results)
  stopCluster(cl)
})



predictedPtype=pRRopheticPredict(
  testMatrix=expdata,
  drug="Erlotinib",
  tissueType="all",
  batchCorrect="eb",
  selection=1,
  dataset="cgp2016")

predictedPtype <- as.data.frame(predictedPtype)
colnames(predictedPtype) <- "Erlotinib"

write.csv(predictedPtype,"Erlotinib.csv")




predictedPtype <- data.table::fread('Erlotinib.csv')

#根据行名新增分组列：
library(stringr)
predictedPtype$Risk <- ifelse(
  str_sub(predictedPtype$V1,1,1)=='L','Low-risk','High-risk'
)
#转换为因子，调整顺序：
predictedPtype$Risk <- factor(predictedPtype$Risk,levels = c('Low-risk','High-risk'))
my_comparisons <- list( c("Low-risk", "High-risk")) #添加比较分组
library(ggplot2)
library(ggpubr)
p1 <- ggviolin(predictedPtype, x = 'Risk', y = 'Erlotinib', fill = 'Risk',
               palette = c("#4979b6","#d9352a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size=0.5, tip.length = 0.02, method = 't.test')
p1
