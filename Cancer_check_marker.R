#############################################2018-2019-merge-Science reference
rm(list = ls())
setwd("/data5/")
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(scater)
library(SingleCellExperiment)
library(dplyr) 
library(ggplot2)
library(gridExtra)
setwd("/data5/zhangq/scRNA/MultiSample/NCTMerge9")
load("Kidney_NCT_merge.rda")
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/C_Anno')
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")

selectobject=ifnb.list[["Cancer"]]
GetAssay(selectobject,assay = "RNA")

######
P0=FeaturePlot(selectobject,features = c("EPCAM", "KRT18", "KRT19","DEFB1","CTSK","SOX4")
               ,reduction="umap"
               ,pt.size =0.4,ncol =3,cols = c("lightgrey","#ff0000"))
ggsave(P0,filename = "Kidney.integrated_Cancer_Epi-check.png",height = 9,width = 15)


P2=FeaturePlot(selectobject,features = c("EPCAM", "PTPRC", "PECAM1","ACTA2")
               ,reduction="umap"
               ,pt.size =0.4,ncol =2,cols = c("lightgrey","#ff0000"))
ggsave(P2,filename = "Kidney.integrated_Cancer_EPCAM-PTPRC-PECAM1-ACTA2.png",height = 16,width = 18)


P22=FeaturePlot(selectobject,features = c("EPCAM","ACTA2"
                                          ,"PECAM1","TPSAB1" #MAST 
                                          ,"CD3D","CD3G" #T
                                          ,"LYZ","C1QB" #Myeloid
                                          ,"CD79A","MS4A1" #B
                                          ,"KLRF1","KLRD1" #NK
                                          ,"IGLL1","MZB1" #Plasma
                                          )
               ,reduction="umap"
               ,pt.size =0.4,ncol =4,cols = c("lightgrey","#ff0000"))
ggsave(P22,filename = "Kidney.integrated_Cancer_new.png",height = 19,width = 22)

#######cancer
Pc=FeaturePlot(selectobject,features = c("NDUFA4L2","CA9","SLC17A3","NNMT")
               ,reduction="umap",ncol=4,
               ,pt.size =0.4,cols = c("lightgrey","#ff0000"))
ggsave(Pc,filename = "Kidney.integrated_Cancer_ccRCC.png",height = 4.5,width = 20)

######Bcell#######
P3=FeaturePlot(selectobject,features = c("CD79A","MS4A1")
               ,reduction="umap"
               ,pt.size =0.4,ncol =2,cols = c("lightgrey","#ff0000"))
ggsave(P3,filename = "Kidney.integrated_Cancer_Bcell.png",height = 4.5,width = 10)
#############T NK
P32=FeaturePlot(selectobject,features = c("CD8A","CD8B","CD3D","CD3G","CD4","IL7R","FOXP3","TOP2A","KLRF1","KLRD1")
                ,reduction="umap"
                ,pt.size =0.4,ncol =4,cols = c("lightgrey","#ff0000"))
ggsave(P32,filename = "Kidney.integrated_Cancer_define_TNK.png",height = 14,width = 22)


pdf(file="Cancer_naiveT_marker.pdf", width = 30, height =7)
pp = DotPlot(selectobject,idents=c("1","3","8","9","12","21")
             ,features = c("PTPRC","CD69","CD28","CD40LG","IL27","IL10"
             ,"RORC","TNFRSF18","TBX21","CX3CR1","GZMB"
             ,"TRGC2","TRDC","KLRB1","SLC4A10","ZBTB16","CCR5"
             ,"CXCR6","EOMES","GZMK","IL2RA","CTLA4","FOXP3","CXCR5"
             ,"CXCR3","SELL","CCR7","CD19","CD14","MKI67","IL7R"
             ,"FCGR3A","KLRF1","CD8B","CD8A","CD4","CD3D"
             ),group.by = 'seurat_clusters'
             #,scale.by="size"
             ,assay = "RNA",cols = c('white',"#ff0000"),dot.scale =12) 
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + RotatedAxis()+
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6),panel.border = element_rect(colour = "black", fill=NA, size=0.6))
print(pp)
dev.off()
ggsave(pp,filename = "Cancer_naiveT_marker.png",height = 6,width =15)


#############Monocyte
P41=FeaturePlot(selectobject,features = c("CD14","FCGR3A","S100A8","S100A9") 
                ,reduction="umap"
                ,pt.size =0.4,ncol =2,cols = c("lightgrey","#ff0000"))
ggsave(P41,filename = "Kidney.integrated_Cancer_define_mono.png",height = 9,width = 10)

##############################macrophage "APOE","APOC1","MRC1"
P42=FeaturePlot(selectobject,features = c("LYZ","CD163","CD68","HLA-DRA") 
                ,reduction="umap"
                ,pt.size =0.4,ncol =2,cols = c("lightgrey","#ff0000"))
ggsave(P42,filename = "Kidney.integrated_Cancer_define_macro.png",height = 9,width = 10)
##############################DC
P43=FeaturePlot(selectobject,features = c("CD1C","CLEC9A") 
                ,reduction="umap"
                ,pt.size =0.4,ncol =2,cols = c("lightgrey","#ff0000"))
ggsave(P43,filename = "Kidney.integrated_Cancer_define_DC.png",height = 4.5,width = 10)

####################################proliferatingT
P5=FeaturePlot(selectobject,features = c("STMN1","MKI67","TOP2A")
               ,reduction="umap"
               ,pt.size =0.4,ncol =3,cols = c("lightgrey","#ff0000"))
ggsave(P5,filename = "Kidney.integrated_Cancer_define_proliferatingT.png",height = 4.5,width = 15)

###########################plasma
P7=FeaturePlot(selectobject,features = c("IGHG3","MZB1","JCHAIN")
               ,reduction="umap"
               ,pt.size =0.4,ncol =3,cols = c("lightgrey","#ff0000"))
ggsave(P7,filename = "Kidney.integrated_Cancer_define_plasma.png",height = 4.5,width = 15)

###########################Neutrophil
P8=FeaturePlot(selectobject,features = c("SELL","CXCL8","FCGR3B","MNDA")
               ,reduction="umap"
               ,pt.size =0.4,ncol =2,cols = c("lightgrey","#ff0000"))
ggsave(P8,filename = "Kidney.integrated_Cancer_define_Neutrophil.png",height = 9,width = 10)

###############################
pp = DotPlot(selectobject,idents=c("0","6","7","10","13","26")
             ,features = c("EPCAM","PTPRC","PECAM1","NDUFA4L2","CA9","SLC17A3","NNMT"),group.by = 'seurat_clusters'
             #,scale.by="size"
             ,assay = "RNA",cols = c('white','#ff0000'),dot.scale =12) 
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6),panel.border = element_rect(colour = "black", fill=NA, size=0.6))
ggsave(pp,filename = "Cancer_ccRCC_check_marker.png",height = 6,width =15)


RCCmarker=list(
        "ccRCC"=c("NDUFA4L2","CA9","SLC17A3","NNMT"),
        "sRCC"=c("PAX8","MUC1","MME","VIM"),
        "EMT"=c("TGFB1","SNAI1","SNAI2","SNAI3"),
        "Myofibroblast"=c("ACTA2","RGS5")
)

pdf(file="Cancer_ccRCC_sRCC_marker.pdf", width = 15, height = 6)
pp = DotPlot(selectobject,idents=c("0","6","7","10","13","26")
             ,features = RCCmarker,group.by = 'seurat_clusters'
             #,scale.by="size"
             ,assay = "RNA",cols = c('white','#ff0000'),dot.scale =12) 
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6),panel.border = element_rect(colour = "black", fill=NA, size=0.6))

print(pp)
dev.off()
ggsave(pp,filename = "Cancer_ccRCC_sRCC_marker.png",height = 6,width =15)
############################macrophage
pp = DotPlot(selectobject,idents=c("2","15","30","31")
             ,features = c("LYZ","CD163","CD68","HLA-DRA","CD14","FCGR3A","S100A8","S100A9","CD1C","CLEC9A"),group.by = 'seurat_clusters'
             #,scale.by="size"
             ,assay = "RNA",cols = c('white',"#ff0000"),dot.scale =12) 
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + RotatedAxis()+
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6),panel.border = element_rect(colour = "black", fill=NA, size=0.6))
ggsave(pp,filename = "Cancer_macrophage_marker.png",height = 6,width =15)

###################################################
pp = DotPlot(selectobject,idents=c("1","3","8","9","12","21")
             ,features = c("PTPRC","CD69","CD28","CD40LG","IL27","IL10"
             ,"RORC","TNFRSF18","TBX21","CX3CR1","GZMB"
             ,"TRGC2","TRDC","KLRB1","SLC4A10","ZBTB16","CCR5"
             ,"CXCR6","EOMES","GZMK","IL2RA","CTLA4","FOXP3","CXCR5"
             ,"CXCR3","SELL","CCR7","CD19","CD14","MKI67","IL7R"
             ,"FCGR3A","KLRF1","CD8B","CD8A","CD4","CD3D"
             ),group.by = 'seurat_clusters'
             #,scale.by="size"
             ,assay = "RNA",cols = c('white','#ff0000'),dot.scale =12) 
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + RotatedAxis()+
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6),panel.border = element_rect(colour = "black", fill=NA, size=0.6))
ggsave(pp,filename = "Cancer_naiveT_marker.png",height = 6,width =15)

################
pp = DotPlot(selectobject,idents=c("8")
             ,features = c("PTPRC","CD69","CD28","CD40LG","IL27","IL10"
             ,"RORC","TNFRSF18","TBX21","CX3CR1","GZMB"
             ,"TRGC2","TRDC","KLRB1","SLC4A10","ZBTB16","CCR5"
             ,"CXCR6","EOMES","GZMK","IL2RA","CTLA4","FOXP3","CXCR5"
             ,"CXCR3","SELL","CCR7","CD19","CD14","MKI67","IL7R"
             ,"FCGR3A","KLRF1","CD8B","CD8A","CD4","CD3D"
             ),group.by = 'seurat_clusters'
             #,scale.by="size"
             ,assay = "RNA",cols = c('white','#ff0000'),dot.scale =12) 
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + RotatedAxis()+
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6),panel.border = element_rect(colour = "black", fill=NA, size=0.6))

ggsave(pp,filename = "Cancer_Cluster_naiveT_marker.png",height = 6,width =15)

###########
###########################Coexpression
P9=FeaturePlot(selectobject,features = c("CD8A","CD4")
               ,reduction="umap"
               ,blend=T
               ,pt.size =0.4)
ggsave(P9,filename = "Cancer_CD8A_CD4_Coexpression.png",height = 4.5,width = 17)

P10=FeaturePlot(selectobject,features = c("PTPRC","CD3D")
               ,reduction="umap"
               ,blend=T
               ,pt.size =0.4)
ggsave(P10,filename = "Cancer_PTPRC_CD3D_Coexpression.png",height = 4.5,width = 17)

P11=FeaturePlot(selectobject,features = c("PTPRC","CD3G")
               ,reduction="umap"
               ,blend=T
               ,pt.size =0.4)
ggsave(P11,filename = "Cancer_PTPRC_CD3G_Coexpression.png",height = 4.5,width = 17)

################## cluster number
cluster_cellnum=as.matrix(table(selectobject@meta.data$seurat_clusters))
cluster_cellnum=cbind(rownames(cluster_cellnum),cluster_cellnum)
colnames(cluster_cellnum)=c("Cluster number","Cell number")
write.csv(cluster_cellnum,"Cancer_cluster_cellnum.csv",row.names = F,quote = F)

#################################################################Endothelial-check
pp = DotPlot(selectobject,idents=c("0","6","7","10","13","26")
             ,features = c("EPCAM","ACTA2","PECAM1","PTPRB","KDR","NDUFA4L2","CA9","SLC17A3","NNMT"),group.by = 'seurat_clusters'
             #,scale.by="size"
             ,assay = "RNA",cols = c('white','#ff0000'),dot.scale =12) 
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6),panel.border = element_rect(colour = "black", fill=NA, size=0.6))
ggsave(pp,filename = "Cancer_Endothelial-check_marker.png",height = 4.5,width =10)


P12=FeaturePlot(selectobject,features = c("PECAM1","EPCAM")
               ,reduction="umap"
               ,blend=T
               ,pt.size =0.4)
ggsave(P12,filename = "Cancer_PECAM1_ACTA2_Coexpression.png",height = 4.5,width = 17)
