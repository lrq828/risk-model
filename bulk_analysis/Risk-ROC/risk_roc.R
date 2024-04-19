library(ggrisk)
library(survival)
library(rms)
data = read.delim("4miRNA-bluk&survial.txt")
rownames(data)=data[,1]
data=data[,-1]

#以hsa.mir.100+hsa.mir.125b.2+hsa.mir.145+hsa.mir.221进行分析
# data=data[,1:6]
fit <- cph(Surv(OS.time,OS)~hsa.mir.100+hsa.mir.125b.2+hsa.mir.145+hsa.mir.221,data)

pdf("风险联动图.pdf")
ggrisk(fit,cutoff.show = F)
dev.off()
write.csv(fit$linear.predictors,file = "TCGA-riskscore.csv")



############## 内部验证########

data = read.delim("4miRNA-bluk&survial.txt")
rownames(data)=data[,1]
data=data[,-1]
set.seed(2023)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.5, 0.5))
data1 <- data[ind == 1,]  # 训练数据
data2 <- data[ind == 2,]  # 测试数据


fit1 <- cph(Surv(OS.time,OS)~hsa.mir.100+hsa.mir.125b.2+hsa.mir.145+hsa.mir.221,data1)

pdf("风险联动图-BIT1.pdf")
ggrisk(fit1,cutoff.show = F)
dev.off()
write.csv(fit1$linear.predictors,file = "BIT-riskscore1.csv")

fit2 <- cph(Surv(OS.time,OS)~hsa.mir.100+hsa.mir.125b.2+hsa.mir.145+hsa.mir.221,data2)

pdf("风险联动图-BIT2.pdf")
ggrisk(fit2,cutoff.show = F)
dev.off()
write.csv(fit2$linear.predictors,file = "BIT-riskscore2.csv")





###################### ROC #####################
library(timeROC);library(caret)
# rm(list = ls())
data <- read.csv("TCGA_Risk.csv",row.names = 1)
dim(data)


train <- data

col <- c("#0073C2FF","firebrick1","orange","yellow","red") ## 自定义颜色
ROC <- timeROC(train$OS.time,train$OS,train$riskScore,
              cause = 1,weighting = 'marginal',
              times=round(quantile(train$OS.time,probs=c(0.42,0.89,0.97)),0),
              ROC = T,iid = T)
(ROC[["times"]])/365
ROC$AUC
confint(ROC)$CI_AUC
pdf("135year(TCGA).pdf")
plot(ROC,time=396,title=FALSE,lwd=1.5,col=col[1])
plot(ROC,time=1130,col=col[2],add=TRUE,title=FALSE,lwd=1.5)
plot(ROC,time=1878,col=col[3],add=TRUE,title=FALSE,lwd=1.5)
id <- c(paste0("Year-1 AUC : ",round(ROC$AUC[1],4),"(0.560-0.701)"),
        paste0("Year-3 AUC : ",round(ROC$AUC[2],4),"(0.543-0.734)"),
        paste0("Year-5 AUC : ",round(ROC$AUC[3],4),"(0.528-0.891)")
)
legend("bottomright",id,
       fill=col[1:3],
       bty="o",cex=1,
       border = NA)
abline(0,1,lty=2,lwd=0.5)
title(main = "TCGA")
dev.off()




############## 内部验证########
rm(list = ls())
data <- read.csv("BIT1_Risk.csv",row.names = 1)
dim(data)

train <- data

#训练集
## ROC
col <- c("#0073C2FF","firebrick1","orange","yellow","red") ## 自定义颜色
ROC <- timeROC(train$OS.time,train$OS,train$riskScore,
               cause = 1,weighting = 'marginal',
               times=round(quantile(train$OS.time,probs=c(0.40,0.9,0.97)),0),
               ROC = T,iid = T)
(ROC[["times"]])/365
ROC$AUC
confint(ROC)$CI_AUC
pdf("135year(TCGA)-BIT1.pdf")
plot(ROC,time=378,title=FALSE,lwd=1.5,col=col[1])
plot(ROC,time=1105,col=col[2],add=TRUE,title=FALSE,lwd=1.5)
plot(ROC,time=1830,col=col[3],add=TRUE,title=FALSE,lwd=1.5)
id <- c(paste0("Year-1 AUC : ",round(ROC$AUC[1],4),"(0.515-0.723)"),
        paste0("Year-3 AUC : ",round(ROC$AUC[2],4),"(0.552-0.798)"),
        paste0("Year-5 AUC : ",round(ROC$AUC[3],4),"(0.6-0.903)")
)
legend("bottomright",id,
       fill=col[1:3],
       bty="o",cex=1,
       border = NA)
abline(0,1,lty=2,lwd=0.5)
title(main = "BIT set1")
dev.off()


######

rm(list = ls())
data <- read.csv("BIT2_Risk.csv",row.names = 1)
dim(data)

train <- data
# validation <- read.delim('')
#训练集
## ROC
col <- c("#0073C2FF","firebrick1","orange","yellow","red") ## 自定义颜色
ROC <- timeROC(train$OS.time,train$OS,train$riskScore,
               cause = 1,weighting = 'marginal',
               times=round(quantile(train$OS.time,probs=c(0.35,0.85,0.96)),0),
               ROC = T,iid = T)
(ROC[["times"]])/365
ROC$AUC
confint(ROC)$CI_AUC
pdf("135year(TCGA)-BIT2.pdf")
plot(ROC,time=389,title=FALSE,lwd=1.5,col=col[1])
plot(ROC,time=1098,col=col[2],add=TRUE,title=FALSE,lwd=1.5)
plot(ROC,time=1882,col=col[3],add=TRUE,title=FALSE,lwd=1.5)
id <- c(paste0("Year-1 AUC : ",round(ROC$AUC[1],4),"(0.546-0.745)"),
        paste0("Year-3 AUC : ",round(ROC$AUC[2],4),"(0.510-0.751)"),
        paste0("Year-5 AUC : ",round(ROC$AUC[3],4),"(0.540-0.940)")
)
legend("bottomright",id,
       fill=col[1:3],
       bty="o",cex=1,
       border = NA)
abline(0,1,lty=2,lwd=0.5)
title(main = "BIT set2")
dev.off()
