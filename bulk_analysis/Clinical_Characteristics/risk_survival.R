library(survival)
library(survminer)
#install.packages("survival")


############## Age(大于等于65)##########
sur <- read.delim('Age(大于等于65).txt',row.names = 1)#临床数据
#使用coxph函数拟合cox风险回归模型，可以使用连续型变量，也可以使用分类型变量（二分类或者多分类都可）
##对年龄（离散型变量）进行cox回归：
cox_Group <- coxph(Surv(OS.time, OS) ~ Group , data = sur)
#查看p值：
cox_Group

summary(cox_Group)
#画生存分析曲线
f <- survfit(Surv(OS.time, OS) ~ Group, data = sur)
summary(f)
p1 <- ggsurvplot(
  f,
  xlab = "Time_days",
  pval = T,
  pval.method = T,
  legend = c(0.85, 0.87),
  title = "Age >= 65",
  legend.title = "Treatment",
  conf.int = T,
  conf.int.style = 'step',
  surv.median.line = 'hv',
  palette = 'lancet',
  risk.table = T,
  risk.table.pos = c('in'),
  ggtheme = theme_bw(base_size = 12)
)


############## Age(小于65)##########
sur <- read.delim('Age(小于65).txt',row.names = 1)#临床数据
#使用coxph函数拟合cox风险回归模型，可以使用连续型变量，也可以使用分类型变量（二分类或者多分类都可）
##对年龄（离散型变量）进行cox回归：
cox_Group <- coxph(Surv(OS.time, OS) ~ Group , data = sur)
#查看p值：
cox_Group

summary(cox_Group)
#画生存分析曲线
f <- survfit(Surv(OS.time, OS) ~ Group, data = sur)
summary(f)
p2 <- ggsurvplot(
  f,
  xlab = "Time_days",
  pval = T,
  pval.method = T,
  legend = c(0.85, 0.87),
  title = "Age < 65",
  legend.title = "Treatment",
  conf.int = T,
  conf.int.style = 'step',
  surv.median.line = 'hv',
  palette = 'lancet',
  risk.table = T,
  risk.table.pos = c('in'),
  ggtheme = theme_bw(base_size = 12)
)


############## Female##########
sur <- read.delim('Female.txt',row.names = 1)#临床数据
#使用coxph函数拟合cox风险回归模型，可以使用连续型变量，也可以使用分类型变量（二分类或者多分类都可）
##对年龄（离散型变量）进行cox回归：
cox_Group <- coxph(Surv(OS.time, OS) ~ Group , data = sur)
#查看p值：
cox_Group

summary(cox_Group)
#画生存分析曲线
f <- survfit(Surv(OS.time, OS) ~ Group, data = sur)
summary(f)
p3 <- ggsurvplot(
  f,
  xlab = "Time_days",
  pval = T,
  pval.method = T,
  legend = c(0.85, 0.87),
  title = "Female",
  legend.title = "Treatment",
  conf.int = T,
  conf.int.style = 'step',
  surv.median.line = 'hv',
  palette = 'lancet',
  risk.table = T,
  risk.table.pos = c('in'),
  ggtheme = theme_bw(base_size = 12)
)


############## Male##########
sur <- read.delim('Male.txt',row.names = 1)#临床数据
#使用coxph函数拟合cox风险回归模型，可以使用连续型变量，也可以使用分类型变量（二分类或者多分类都可）
##对年龄（离散型变量）进行cox回归：
cox_Group <- coxph(Surv(OS.time, OS) ~ Group , data = sur)
#查看p值：
cox_Group

summary(cox_Group)
#画生存分析曲线
f <- survfit(Surv(OS.time, OS) ~ Group, data = sur)
summary(f)
p4 <- ggsurvplot(
  f,
  xlab = "Time_days",
  pval = T,
  pval.method = T,
  legend = c(0.85, 0.87),
  title = "Male",
  legend.title = "Treatment",
  conf.int = T,
  conf.int.style = 'step',
  surv.median.line = 'hv',
  palette = 'lancet',
  risk.table = T,
  risk.table.pos = c('in'),
  ggtheme = theme_bw(base_size = 12)
)


############## Stage(Ⅰ-Ⅱ)##########
sur <- read.delim('Stage(Ⅰ-Ⅱ).txt',row.names = 1)#临床数据
#使用coxph函数拟合cox风险回归模型，可以使用连续型变量，也可以使用分类型变量（二分类或者多分类都可）
##对年龄（离散型变量）进行cox回归：
cox_Group <- coxph(Surv(OS.time, OS) ~ Group , data = sur)
#查看p值：
cox_Group

summary(cox_Group)
#画生存分析曲线
f <- survfit(Surv(OS.time, OS) ~ Group, data = sur)
summary(f)
p5 <- ggsurvplot(
  f,
  xlab = "Time_days",
  pval = T,
  pval.method = T,
  legend = c(0.85, 0.87),
  title = "Stage(Ⅰ-Ⅱ)",
  legend.title = "Treatment",
  conf.int = T,
  conf.int.style = 'step',
  surv.median.line = 'hv',
  palette = 'lancet',
  risk.table = T,
  risk.table.pos = c('in'),
  ggtheme = theme_bw(base_size = 12)
)


############## Stage(Ⅲ-Ⅳ).#########
sur <- read.delim('Stage(Ⅲ-Ⅳ).txt',row.names = 1)#临床数据
#使用coxph函数拟合cox风险回归模型，可以使用连续型变量，也可以使用分类型变量（二分类或者多分类都可）
##对年龄（离散型变量）进行cox回归：
cox_Group <- coxph(Surv(OS.time, OS) ~ Group , data = sur)
#查看p值：
cox_Group

summary(cox_Group)
#画生存分析曲线
f <- survfit(Surv(OS.time, OS) ~ Group, data = sur)
summary(f)
p6 <- ggsurvplot(
  f,
  xlab = "Time_days",
  pval = T,
  pval.method = T,
  legend = c(0.85, 0.87),
  title = "Stage(Ⅲ-Ⅳ)." , 
  legend.title = "Treatment",
  conf.int = T,
  conf.int.style = 'step',
  surv.median.line = 'hv',
  palette = 'lancet',
  risk.table = T,
  risk.table.pos = c('in'),
  ggtheme = theme_bw(base_size = 12)
)



pdf("临床特征生存分析曲线.pdf")
p1;p2;p3;p4;p5;p6
dev.off()
library(patchwork)
# ggsave(p,filename = "临床特征生存分析曲线.pdf",width = 15,height =12)
par(mfrow=c(3,2))
p1;p2;p3;p4;p5;p6



library(ggpubr)

ggarrange(p1,p2,p3,p4,ncol=2,nrow=2,labels=c("A","B","C","D"))
saveRDS(c(p1,p2,p3,p4,p5,p6),"plot.RDS")
