setwd("E:/GC/bluk/免疫治疗分析/TIDE&EaSIeR/")

############ FPKM转TPM ##########
# expMatrix <- read.table(file="genes_exp_fpkm(331samples).txt",header=TRUE, sep="\t",check.names=F,quote="",row.names = 1)
# 
# fpkmToTpm <- function(fpkm)
# {
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# tpms_expMatrix <- apply(expMatrix,2,fpkmToTpm)
# tpms_expMatrix[1:3,]
# colSums(tpms_expMatrix)
# 
# write.csv(tpms_expMatrix,"tpms_expMatrix(331samples).csv")

###########


library("easier")
RNA_tpm <- read.csv('tpms_expMatrix(331samples).csv',row.names = 1)
{
  # 考虑到部分基因有多个对应关系，需要进一步处理（保留作者给定的gene symbol）
  genes_info <- easier:::reannotate_genes(cur_genes = rownames(RNA_tpm))
  
  ## 去除不支持的基因symbol
  non_na <- !is.na(genes_info$new_names)
  RNA_tpm <- RNA_tpm[non_na, ]
  genes_info <- genes_info[non_na, ]
  
  ## 去除 entries that are withdrawn
  RNA_tpm <- RNA_tpm[-which(genes_info$new_names == "entry withdrawn"), ]
  genes_info <- genes_info[-which(genes_info$new_names == "entry withdrawn"), ]
  
  ## 找出重复基因
  # newnames_dup <- unique(genes_info$new_names[duplicated(genes_info$new_names)])
  # newnames_dup_ind <- do.call(c, lapply(newnames_dup, function(X) which(genes_info$new_names == X)))
  # newnames_dup <- genes_info$new_names[newnames_dup_ind]
  # 
  # ## 检索重复基因的数据
  # tmp <- RNA_tpm[genes_info$old_names[genes_info$new_names %in% newnames_dup],]
  
  ## 删除重复基因的数据
  # RNA_tpm <- RNA_tpm[-which(rownames(RNA_tpm) %in% rownames(tmp)),]
  # 
  # ## 整合重复基因的数据
  # dup_genes <- genes_info$new_names[which(genes_info$new_names %in% newnames_dup)]
  # names(dup_genes) <- rownames(tmp)
  # if (anyDuplicated(newnames_dup)){
  #   tmp2 <- stats::aggregate(tmp, by = list(dup_genes), FUN = "mean")
  #   rownames(tmp2) <- tmp2$Group.1
  #   tmp2$Group.1 <- NULL
  # }
  # 
  # # 整理归纳到一个表达矩阵
  # RNA_tpm <- rbind(RNA_tpm, tmp2)
  

}

hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
immune_response_scores <- compute_scores_immune_response(RNA_tpm = RNA_tpm, 
                                                         selected_scores = hallmarks_of_immune_response)
head(immune_response_scores) 
###########细胞成分######
cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)


###########通路活性方面######
RNA_counts <- read.csv('counts_expMatrix(331samples).csv',row.names = 1)
pathway_activities <- compute_pathway_activity(RNA_counts = RNA_counts,
                                               remove_sig_genes_immune_response = TRUE)
head(pathway_activities)


###########转录因子活性######
tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)
head(tf_activities[,1:5])


###########配体-受体权重与细胞间相互作用######
lrpair_weights <- compute_LR_pairs(RNA_tpm = RNA_tpm,
                                   cancer_type = "pancan")
head(lrpair_weights)[,1:5]

ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, 
                                  cancer_type = "pancan")
head(ccpair_scores) 
# names(intercell_networks)




###############预测患者的免疫反应###########
predictions <- predict_immune_response(pathways = pathway_activities,
                                       immunecells = cell_fractions,
                                       tfs = tf_activities,
                                       lrpairs = lrpair_weights,
                                       ccpairs = ccpair_scores,
                                       cancer_type = "STAD", 
                                       verbose = TRUE)
# predict_immune_response


available_cancer_types <- c("STAD")


TMB.all <- read.delim("TMB.txt",row.names = 1)
TMB <- TMB.all$total
names(TMB) <- row.names(TMB.all)
################### 获得免疫反应评分########################
easier_derived_scores <- retrieve_easier_score(predictions_immune_response = predictions,
                                               TMB_values = TMB,
                                               easier_with_TMB = c("weighted_average", 
                                                                   "penalized_score"),
                                               weight_penalty = 0.5)

save(cell_fractions,easier_derived_scores,file = "easier_derived_scores&cell-fraction.Rdata")




###################使用有免疫治疗反应数据评估EaSIeR的预测########################
output_eval_with_resp <- assess_immune_response(predictions_immune_response = predictions,
                                                patient_response = patient_ICBresponse,
                                                RNA_tpm = RNA_tpm,
                                                TMB_values = TMB,
                                                easier_with_TMB = "weighted_average",
                                                weight_penalty = 0.5)


###################无免疫治疗数据评估免疫反应########################
output_eval_no_resp <- assess_immune_response(predictions_immune_response = predictions,
                                              TMB_values = TMB,
                                              easier_with_TMB = "weighted_average",
                                              weight_penalty = 0.5)

output_eval_no_resp[[1]]
output_eval_no_resp[[2]]


























