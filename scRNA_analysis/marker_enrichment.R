### 基因ID转换
library(rvcheck)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)


rm(list = ls())
alldegs <- read.csv('memeory B_cell.csv')
# topgenes <- row.names(alldegs[order(alldegs$p_val_adj,decreasing = T)[1:100],])
topgenes <- as.data.frame(alldegs);topgenes <- c(topgenes$genes)
topgenes.fd <- bitr(topgenes,fromType = "SYMBOL",toType = c("ENTREZID","ENSEMBL"),org.Hs.eg.db)
write.csv(topgenes.fd$ENSEMBL,"memeory B_cell_gene_ENSEMBL.csv")

enrich2plot <- function(data){
  library(ggplot2)
  data <- data[order(data$qvalue,decreasing = F)[1:20],]
  data$BgRatio <- 
    apply(data,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])})/apply(data,1,function(x){
        as.numeric(strsplit(x[4],'/')[[1]][2])
      })
  
  p <- ggplot(data,aes(BgRatio,Description))
  p <- p+geom_point()
  
  pbubble <- p + geom_point(aes(size = Count, color = -1*log10(qvalue)))
  
  pr <- pbubble + scale_color_gradient(low = "lightgreen",high = "red")+
    labs(color = expression( -log[10](qvalue)),size = "observer.gene.count",x = "Richfactor",
         y = "term.description",title = "Enrichment Result")
  
  pr <- pr + theme_bw()
  pr
}

########kegg
kegg <- enrichKEGG(topgenes.fd$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
kegg <- setReadable(kegg,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID')
write.csv(kegg@result,'B_memory_kegg_result.csv')

pdf('B_memory_kegg.pdf',width = 10,height = 10)
dotplot(kegg)
dev.off()

pdf('B_memory_kegg_barplot.pdf',width = 10,height = 10)
barplot(kegg)
dev.off()

########GO
go <- enrichGO(topgenes.fd$SYMBOL, org.Hs.eg.db,keyType = "SYMBOL", ont = "ALL",
                 pvalueCutoff = 0.05,qvalueCutoff = 0.05)

write.csv(go@result,'B_memory_go_result.csv')
pdf('B_memory_go.pdf',width = 10,height = 10)
dotplot(go)
dev.off()

go_gene <- go@gene
go_gene <- bitr(go_gene,fromType = "SYMBOL",toType = c("ENSEMBL"),org.Hs.eg.db)
write.csv(go_gene,"go_gene.csv")
kegg_gene <- kegg@gene
kegg_gene <- bitr(kegg_gene,fromType = "ENTREZID",toType = c("ENSEMBL"),org.Hs.eg.db)
write.csv(kegg_gene,"kegg_gene.csv")
