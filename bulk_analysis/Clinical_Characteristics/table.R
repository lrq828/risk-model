# install.packages("compareGroups")  # 安装包
library(compareGroups)  # 加载包
data <- read.delim("age_gender_stage_TNM.txt",row.names = 1)
str(data) # 查看数据集结构
descrTable( ~ ., data = data) # 描述下总样本人群

#######选择分组变量
tb <- descrTable(Risk ~ ., data = data)

export2pdf(tb, file='tb.pdf')
