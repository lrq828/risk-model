STADdata=read.table("TCGA-STAD.htseq_counts.tsv",header=T,sep='\t')
STADdata[1:4,1:4]

# 4、去除Ensembl_ID小数点后面的数字
library(tidyr)
STADdata1<-separate(STADdata,Ensembl_ID,into= c("Ensembl_ID"),sep="\\.")
STADdata1[1:4,1:4]

# 5、下载Homo_sapiens.GRCh38.99.chr.gtf.gz文件，然后放到我们R语言的工作目录内，
# 打开文件（可以在网上搜到下载方法，此不在详述），进入相应的工作路径
# setwd("C:\\Users\\86188\\Desktop\\CRC-TCGA\\UCSC xena\\2. Ensemble_ID转换为Gene_Symbol")
# gtf <- rtracklayer::import('Homo_sapiens.GRCh38.109.chr.gtf.gz')
# class(gtf)
# 
# # 6、转换为数据框
# gtf <-as.data.frame(gtf)
# 
# # 7、查看文件，保存文件为Rdata，将来方便我们直接打开
# dim(gtf) 
# save(gtf,file ="Homo_sapiens.GRCh38.109基因组注释文件.Rda")

# 8、由于GTF文件中，基因ID的列名是gene_id,所以我们把STADdata1中的基因列名改成一样的，方便后续合并
colnames(STADdata1)[1]<-'gene_id'

# 9、把蛋白编码的基因的行都筛选出来
a=dplyr::filter(gtf,type=="gene",gene_biotype=="protein_coding")
dim(a)

# 10、只选择gene_name，gene_id和gene_biotype这三列，其他都不要了
b=dplyr::select(a,c(gene_name,gene_id,gene_biotype))
b[1:4,]

# 11、接下来用STADdata1和b文件中共有的gene_id列合并文件
c=dplyr::inner_join(b,STADdata1,by ="gene_id")
c[1:5,1:5]

# 12、再去掉第2列(gene_id)，3列(gene_biotype)，基因名再去重（distinct函数）
library(dplyr)
d=select(c,-gene_id,-gene_biotype)
mRNAdata=distinct(d,gene_name,.keep_all = T)
dim(mRNAdata)

mRNAdata <- mRNAdata[!is.na(mRNAdata[,1]),]                   ## 提取第一列不为空的数据
dim(mRNAdata)


# 13、把行名由数字换成基因
rownames(mRNAdata)<-mRNAdata[,1]
mRNAdata<-mRNAdata[,-1]

# 14、我们下载的数据取过了log2（count+1），这里我们再返回count
mRNAdata<-2^mRNAdata -1
View(mRNAdata)
class(mRNAdata)#此时是数据框“data.frame”，与“matrix”区别就是前者行名不允许有重复的


# 15、本例使用相对更苛刻的原则：若某基因在大于80%的样本中read counts数<10，则将其滤除。
qualified_genes <- c()
for (genes_in_sheet in rownames(mRNAdata)) {
  qualification <- mRNAdata[genes_in_sheet,] <= 10
  if (sum(qualification) < 0.8*length(mRNAdata)) {
    qualified_genes <- append(qualified_genes,genes_in_sheet)
  }
}
mRNAdata <- mRNAdata[qualified_genes,]
dim(mRNAdata)

# 16、保存文件，常见空白分隔符有：sep=” ”；sep = “\t”；sep = “\n”，即空格，制表符，换行符 。到此为止，差异分析前的基因表达数据就整理好了
write.csv(mRNAdata,"mRNAdata_count.csv",quote = F,row.names = T)
# write.table(mRNAdata, file = "mRNAdata.txt", sep = "\t",row.names = T,col.names = NA,quote = F)#col.names = NA表示第1个列名为空
# save(mRNAdata,file ="mRNAdata.Rda")
