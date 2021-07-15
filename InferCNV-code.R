rm(list=ls())
library(infercnv)
library(dplyr)
setwd("/your_path")
load("inferCNV_exprMatrix.rda")
##Create inferCNV object
infercnv_obj = CreateInfercnvObject(delim = '\t',
                                    raw_counts_matrix = exprMatrix,
                                    annotations_file  = '0804_inferCNV_cellAnnota.txt',
                                    gene_order_file   = 'ref_GRCh38_10X_cellranger_GTF_gen_pos_for_inferCNV_gene_name.txt',
                                    ref_group_names = c("Endothelial_cells"
                                                        ,"Myofibroblasts"))

#perform infercnv operations to reveal cnv signal

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='./output_dir',
                             cluster_by_groups=FALSE,
                             k_obs_groups=6,
                             denoise=TRUE,
                             HMM=TRUE,
                             analysis_mode='subclusters')




