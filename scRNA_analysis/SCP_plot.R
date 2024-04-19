library(SCP)

GC <- readRDS("./GC_data/GC_sub_B.rds")
mir_target <- read.table("./GC_data/target.txt",header = T)


GC <- AnnotateFeatures(GC, species = "Homo_sapiens", db = c("TF", "CSPA"))

GC$celltype <- GC@active.ident

ht2 <- FeatureHeatmap(
  srt = GC, group.by = "celltype", features = mir_target$Gene, feature_split = mir_target$Group,
  species = "Homo_sapiens", db = c("GO_BP", "GO_MF", "KEGG"), anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4,topTerm = 8)

pdf("target_FeatureHeatmap3.pdf", height = 10, width = 23)
print(ht2$plot)
dev.off()
