if(!require('multtest'))install.packages("multtest")
if(!require('Seurat'))install.packages("Seurat")
if(!require('dplyr'))install.packages("dplyr")
if(!require('patchwork'))install.packages("patchwork")
if(!require('R.utils'))install.packages("R.utils")
library(ggplot2)
library(gridExtra)


rm(list = ls())
GSM5004180_data <- Read10X(data.dir = "GSE163558/GSM5004180/")
GSM5004181_data <- Read10X(data.dir = "GSE163558/GSM5004181/")
GSM5004182_data <- Read10X(data.dir = "GSE163558/GSM5004182/")

GSM5004180_s <- CreateSeuratObject(count = GSM5004180_data, project = "GSM5004180", min.cells=3, min.features=500)
GSM5004181_s <- CreateSeuratObject(count = GSM5004181_data, project = "GSM5004181", min.cells=3, min.features=500)
GSM5004182_s <- CreateSeuratObject(count = GSM5004182_data, project = "GSM5004182", min.cells=3, min.features=500)

GC<- merge(GSM5004180_s, y = c(GSM5004181_s,GSM5004182_s),
                     add.cell.ids = c("S1", "S2","S3"), project = "GSE163558")

GC
dim(GC)
GC[["percent.mt"]] <- PercentageFeatureSet(GC, pattern = "^MT-")

# pdf('A.pdf',width = 15,height = 10)
VlnPlot(GC, features = c("nFeature_RNA","nCount_RNA", "percent.mt"),ncol=3)
# dev.off()


GC <- subset(GC, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

pdf('QC.pdf',width = 15,height = 10)
VlnPlot(GC, features = c("nFeature_RNA","nCount_RNA", "percent.mt"),ncol=3)
dev.off()

GC <- NormalizeData(GC, normalization.method = "LogNormalize", scale.factor = 10000)
GC <- FindVariableFeatures(GC, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GC),10)
plot1 <- VariableFeaturePlot(GC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf('VariableFeature.pdf',width = 10,height = 10)
plot2 
dev.off()

GC <- ScaleData(GC, features = rownames(GC))
GC <- RunPCA(GC, features = VariableFeatures(object = GC))

VizDimLoadings(GC,dims = 1:5, reduction = "pca")
#投影的降维图
DimPlot(GC, reduction = "pca")
# 前15个pca维度的表现
DimHeatmap(GC, dims = 1:20,cells = 500,balanced = TRUE) # 方法一

### 高阶pca分析
GC <- JackStraw(GC, num.replicate = 100)
GC <- ScoreJackStraw(GC, dims = 1:20)
pdf('JackStraw.pdf',width = 15,height = 10)
JackStrawPlot(GC, dims = 1:20) # 方法二
dev.off()

#考虑选择多少个主成分，鉴于速度问题，性价比高的方式是用碎石图
pdf('Elbow.pdf',width = 10,height = 10)
ElbowPlot(GC) 
dev.off()



GC <- FindNeighbors(GC, dims = 1:20)
GC <- FindClusters(GC, resolution = 0.1)
GC <- RunTSNE(GC, dims = 1:20)
DimPlot(GC, reduction = "tsne")
GC <- RunUMAP(GC, dims = 1:20)
# 各样本的情况
# DimPlot(GC, reduction = "umap", group.by = "orig.ident")
pdf('umap_unannotation.pdf',width = 10,height = 10)
DimPlot(GC, reduction = "umap")
dev.off()




#######查看各群maker表达#########
# a <- list(Monocyte=c("S100A8","S100A9","CD14"),
#           B=c("MS4A1","LINC00926","CD79A"),
#           Plasma=c("MZB1","IGKC","JCHAIN"),
#           NK=c("GNLY","NKG7","KLRD1"),
#           T=c("CD3D","CD3E"),
#           DC=c("IL3RA","GZMB","CST3"),
#           Epilthelium=c("KRT5","KRT14"),
#           Myeloid=c("FCGR3A","ITGAM")
#           )
b <- list(B=c("MS4A1","IGHG1","CD79A"),
          Epilthelium=c("MUC1","CLDN4","KRT18","KRT19","EPCAM"),
          Stromal=c("MCAM","VWF","DCN","COL1A2"),
          Proliferative=c("PCNA","MKI67","STMN1"),
          Myeloid=c("CD68","CSF1R","CSF3R"),
          T=c("CD3D","CD3E","CD3G","CD2"),
          NK=c("GNLY","NKG7","KLRD1"))

pdf('DotPlot_Maker.pdf',width = 20,height = 10)
DotPlot(GC,features=b)+RotatedAxis()+labs(x='',y='')
dev.off()


#######结合DotPlot结果去除11群######
GC <- GC1
GC <- subset(GC, idents = c("0","1","2","3","4","5","6","7","8","9","10"))

pdf('DotPlot_Maker.pdf',width = 20,height = 10)
DotPlot(GC,features=b)+RotatedAxis()+labs(x='',y='')
dev.off()

pdf('tsne_unannotation.pdf',width = 10,height = 10)
DimPlot(GC, reduction = "tsne")
dev.off()

pdf('umap_unannotation.pdf',width = 10,height = 10)
DimPlot(GC, reduction = "umap")
dev.off()


###########细胞注释#####
new.cluster.ids <- c("NK cell", "T cell", "Epilthelium cell", "Myeloid cell", "Myeloid cell", "B cell",
                     "Stromal cell", "Proliferative cell", "Stromal cell", "Epilthelium cell", "B cell")
# 把GC的levels传输给new.cluster.ids作为他的name属性
names(new.cluster.ids) <- levels(GC)

GC <- RenameIdents(GC, new.cluster.ids)
levels(GC)

pdf('tsne_annotation.pdf',width = 10,height = 10)
DimPlot(GC, reduction = "tsne", label=T)
dev.off()

pdf('umap_annotation_label.pdf',width = 10,height = 10)
DimPlot(GC, reduction = "umap", label=T)+NoLegend()
dev.off()

b1 <- list(Proliferative=c("PCNA","MKI67","STMN1"),
           Stromal=c("MCAM","VWF","DCN","COL1A2"),
           B=c("MS4A1","IGHG1","CD79A"),
           Myeloid=c("CD68","CSF1R","CSF3R"),
           Epilthelium=c("MUC1","CLDN4","KRT18","KRT19","EPCAM"),
           T=c("CD3D","CD3E","CD3G","CD2"),
           NK=c("GNLY","NKG7","KLRD1"))

pdf('DotPlot_annotation.pdf',width = 20,height = 10)
DotPlot(GC,features=b1)+RotatedAxis()+labs(x='',y='')
dev.off()

saveRDS(GC, file = "./GC.rds")



GC <- readRDS("GC.rds")



#################提取B cell###############
GC_B <- subset(GC, idents = c("B cell"))

GC_B <- FindVariableFeatures(GC_B, selection.method = "vst", nfeatures = 2000)
GC_B <- ScaleData(GC_B)
GC_B <- RunPCA(GC_B, features = VariableFeatures(GC_B))
GC_B <- FindNeighbors(GC_B, dims = 1:20)
GC_B <- FindClusters(GC_B, resolution = 0.4)
GC_B <- RunUMAP(GC_B, dims = 1:20)

DimPlot(GC_B, reduction = "umap")

DotPlot(GC_B,features=c("CD3D","IGHD","IGHG1","LMO2","BCL6","CD27"))

new.cluster.ids <- c("Switched memory B cell", "Naive B cell", "Plasma cell",
                     "Germinal center B cell", "T cell-like B cell")

names(new.cluster.ids) <- levels(GC_B)

GC_B <- RenameIdents(GC_B, new.cluster.ids)
levels(GC_B)

pdf('B_Cell_umap_annotation_label.pdf',width = 10,height = 10)
DimPlot(GC_B, reduction = "umap", label=T)+NoLegend()
dev.off()

pdf('B_Cell_DotPlot_annotation.pdf',width = 10,height = 10)
DotPlot(GC_B,features=c("CD3D","CD27","IGHG1","IGHD","BCL6","LMO2"))+RotatedAxis()+labs(x='',y='')
dev.off()

saveRDS(GC_B, file = "./GC_B.rds")


#######################放回原来大群中########

Idents(GC, cells = colnames(GC_B))<- Idents(GC_B)

GC_sub_B <- GC

pdf('umap_annotation(sub_B).pdf',width = 10,height = 10)
DimPlot(GC_sub_B,label = T)+NoLegend()
dev.off()



saveRDS(GC_sub_B, file = "./GC_sub_B.rds")



GC_sub_B <- readRDS("GC_sub_B.rds")


MB_marker <- FindMarkers(GC_sub_B,ident.1 = "Switched memory B cell",min.pct = 0.3)
write.csv(MB_marker,"MB_Marker.csv")
pdf("net_gene1.pdf", width = 10,height = 10)
VlnPlot(GC_sub_B,c("SETBP1","ZBTB20","TCF4","PKIG","CKS2","HNRNPA2B1"))
dev.off()



#####################UMAP美化#############
{
library(plot1cell)
###Check and see the meta data info on your Seurat object
colnames(GC_sub_B@meta.data)  

###Prepare data for ploting 准备圈图数据
circ_data <- prepare_circlize_data(GC_sub_B, scale = 0.8 )
set.seed(1234)
# 设置细胞分群信息的颜色
cluster_colors<-rand_color(length(levels(GC_sub_B)))
# group_colors<-rand_color(length(names(table(GC_sub_B$orig.ident))))
rep_colors<-rand_color(length(names(table(GC_sub_B$orig.ident))))

###plot and save figures
# 绘制细胞分群圈图
png(filename =  'circlize_plot.png', width = 6, height = 6,units = 'in', res = 300)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
# 添加细胞群注释信息
# add_track(circ_data, group = "Group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()



# 绘制细胞分群圈图
pdf('circlize_plot.pdf', width = 10, height = 10)
plot_circlize(circ_data,do.label = T, pt.size = 0.01, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
# 添加细胞群注释信息
# add_track(circ_data, group = "Group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
dev.off()



library(RColorBrewer)

# 选择调色板
palette <- brewer.pal(length(levels(GC_sub_B)), "Dark2")

# 或者手动设置颜色
# palette <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")

# 生成细胞分群颜色
cluster_colors <- sample(palette, length(levels(GC_sub_B)), replace = TRUE)

# 生成样本颜色
rep_colors <- sample(palette, length(unique(GC_sub_B$orig.ident)), replace = TRUE)

# 保存为 PDF
pdf('circlize_plot2.pdf', width = 6, height = 6)
plot_circlize(circ_data, do.label = TRUE, pt.size = 0.01, col.use = cluster_colors, bg.color = 'white', kde2d.n = 200, repel = TRUE, label.cex = 0.6)
add_track(circ_data, group = "orig.ident", colors = rep_colors, track_num = 3)
dev.off()
}

























