rm(list = ls())
est_dd=read.delim("24miRNA-bluk&survial.txt",row.names = 1)
library(randomForestSRC)
library(survival)

#######RSF#####
rf_nodesize <- 15
seed <- 2022
cindex <- function(data,RS){
  rs=as.data.frame(cbind(data[,1:2],RS))
  aa=coxph(Surv(OS.time,OS)~RS,rs)
  cindex=as.numeric(summary(aa)$concordance[1])
  return(cindex)
}
RSF_fit <- rfsrc(Surv(OS.time,OS)~.,data = est_dd,
             ntree = 1500,nodesize = rf_nodesize,##该值建议多调整  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)

RSF_train_cindex=cindex(est_dd,RSF_fit$predicted)
RSF_miRNA <- names(RSF_fit$importance[RSF_fit$importance>0])
RSF_miRNA <- as.data.frame(RSF_miRNA)
names(RSF_miRNA) <- "miRNA"
pdf("RSF.pdf")
plot(RSF_fit)
dev.off()
write.csv(RSF_miRNA,"RSF-result.csv")

#####StepCox####
StepCox_fit <- step(coxph(Surv(OS.time,OS)~.,est_dd),direction = "both")

StepCox_train_RS = predict(StepCox_fit, type = 'risk', newdata = est_dd)
StepCox_train_cindex=cindex(est_dd,StepCox_train_RS)
StepCox_miRNA  <- names(coef(StepCox_fit))
StepCox_miRNA  <- as.data.frame(StepCox_miRNA) 
names(StepCox_miRNA) <- "miRNA"

write.csv(StepCox_miRNA,"StepCox-result.csv")


# install.packages("ggVennDiagram")
library(ggVennDiagram)
venn_list <- list(RSF_miRNA$miRNA,StepCox_miRNA$miRNA)
pdf("venn.pdf")
ggVennDiagram(venn_list,edge_size=1,set_size = 6,size=5,label_alpha=0,label_size=4,
              category.names = c("RSF","StepCox"))
dev.off()


library(dplyr)
venn <- intersect(RSF_miRNA,StepCox_miRNA)
write.csv(venn,"RSF&StepCox-result.csv")



