#############################################2018-2019-merge-Science reference
rm(list = ls())
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(scater)
library(SingleCellExperiment)
library(dplyr) 
library(ggplot2)
library(gridExtra)
library(cowplot)
setwd("/data5/zhangq/scRNA/thrombus/myRef2")
load("Normal_Merge_Reference.rda")
normal=Project.se
setwd("/data5/zhangq/scRNA/MultiSample/NCTMerge9")
load("Kidney_NCT_merge.rda")
dir.create('/data5/zhangq/scRNA/MultiSample/NCTMerge9/N_Anno')
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/N_Anno')
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")

selectobject=ifnb.list[["Normal"]]
counts <- selectobject@assays$RNA@data 

temp1=as.data.frame(selectobject@active.ident)
temp2=cbind(rownames(temp1),temp1)

colData <- DataFrame(Clusterlabel=temp2[,2],
                     row.names=temp2[,1])

metadata <- "Single cell SummarizedExperiment logCounts" 
Kidney.integrated.se <- SummarizedExperiment(assays=list(logcounts=counts),
                                             colData=colData,metadata=metadata)

set.seed(123)
pred.combined <- SingleR(test = Kidney.integrated.se,
                         ref = normal,
                         labels = normal$label.fine,
                         #de.method="wilcox",
                         assay.type.test = "logcounts",
                         assay.type.ref = "logcounts",
                         method="cluster", clusters=Kidney.integrated.se$Clusterlabel)

table(predicted=pred.combined$labels, truth=rownames(pred.combined))

all.markers <- metadata(pred.combined)$de.genes
pred_mat=as.data.frame(cbind(rownames(pred.combined),pred.combined$labels))
colnames(pred_mat)=c("key","Annotation")
rawcluster=as.data.frame(Kidney.integrated.se$Clusterlabel)
colnames(rawcluster)="key"
newcluster=left_join(rawcluster,pred_mat,by="key")
Kidney.integrated.se$labels <- newcluster[,2]

library(scater)
library(SingleCellExperiment)
Kidney.integrated.sce=as(Kidney.integrated.se,"SingleCellExperiment")

save(pred.combined,Kidney.integrated.sce,Kidney.integrated.se,file = "Kidney.integrated_Normal_pred.rda")


rm(list=ls())
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
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/N_Anno')
load("Kidney.integrated_Normal_pred.rda")
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")
selectobject=ifnb.list[["Normal"]]

pred_mat=as.data.frame(cbind(rownames(pred.combined),pred.combined$labels))
colnames(pred_mat)=c("key","Annotation")

new.cluster.ids <- pred_mat[,2]
names(new.cluster.ids) <- levels(selectobject)

selectobject <- RenameIdents(selectobject, new.cluster.ids)


plot1<-DimPlot(selectobject, reduction = "umap", label = TRUE, pt.size = 0.4) + NoLegend()
plot2<-DimPlot(selectobject, reduction = "tsne", label = TRUE, pt.size = 0.4) + NoLegend()

p4=plot_grid(plot1,plot2
          , labels = c("Science ref", "Science ref")
          , ncol = 2)
ggsave(p4,filename = "Kidney.integrated_Normal_Anno_cellname_umap_tsne.png",height = 6,width = 14)

write.csv(pred_mat,"2018_2019Mergeref_Anno.csv",row.names = F,quote = F)

###############################################
##############################################BP
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
setwd("/data5/zhangq/scRNA/thrombus/loose/")
load("BlueprintEncodeData.rda")
setwd("/data5/zhangq/scRNA/MultiSample/NCTMerge9")
load("Kidney_NCT_merge.rda")
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/N_Anno')
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")
selectobject=ifnb.list[["Normal"]]
counts <- selectobject@assays$RNA@data 

temp1=as.data.frame(selectobject@active.ident)
temp2=cbind(rownames(temp1),temp1)

colData <- DataFrame(Clusterlabel=temp2[,2],
                     row.names=temp2[,1])

metadata <- "Single cell SummarizedExperiment logCounts" 
Kidney.integrated.se <- SummarizedExperiment(assays=list(logcounts=counts),
                                             colData=colData,metadata=metadata)

set.seed(123)
pred.combined <- SingleR(test = Kidney.integrated.se,
                         ref = bp.se,
                         labels = bp.se$label.fine,
                         #de.method="wilcox",
                         assay.type.test = "logcounts",
                         assay.type.ref = "logcounts",
                         method="cluster", clusters=Kidney.integrated.se$Clusterlabel)


table(predicted=pred.combined$labels, truth=rownames(pred.combined))

all.markers <- metadata(pred.combined)$de.genes
pred_mat=as.data.frame(cbind(rownames(pred.combined),pred.combined$labels))
colnames(pred_mat)=c("key","Annotation")
rawcluster=as.data.frame(Kidney.integrated.se$Clusterlabel)
colnames(rawcluster)="key"

newcluster=left_join(rawcluster,pred_mat,by="key")

Kidney.integrated.se$labels <- newcluster[,2]

library(scater)
library(SingleCellExperiment)
Kidney.integrated.sce=as(Kidney.integrated.se,"SingleCellExperiment")


save(pred.combined,Kidney.integrated.sce,Kidney.integrated.se,file = "Kidney.integrated_BP_Normal_pred.rda")

########################################
rm(list=ls())
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(scater)
library(SingleCellExperiment)
library(dplyr) 
library(ggplot2)
library(gridExtra)
library(cowplot)
setwd("/data5/zhangq/scRNA/MultiSample/NCTMerge9")
load("Kidney_NCT_merge.rda")
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/N_Anno')
load("Kidney.integrated_BP_Normal_pred.rda")
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")
selectobject=ifnb.list[["Normal"]]

pred_mat=as.data.frame(cbind(rownames(pred.combined),pred.combined$labels))
colnames(pred_mat)=c("key","Annotation")

new.cluster.ids <- pred_mat[,2]
names(new.cluster.ids) <- levels(selectobject)

selectobject <- RenameIdents(selectobject, new.cluster.ids)
plot1<-DimPlot(selectobject, reduction = "umap", label = TRUE, pt.size = 0.4) + NoLegend()
plot2<-DimPlot(selectobject, reduction = "tsne", label = TRUE, pt.size = 0.4) + NoLegend()

p4=plot_grid(plot1,plot2,ncol = 2
             ,labels = c("SingleR ref","SingleR ref"))
ggsave(p4,filename = "Kidney.integrated_BP_Normal_Anno_cellname_umap_tsne.png",height = 6,width = 14)

write.csv(pred_mat,"BlueprintEncodeData_Anno.csv",row.names = F,quote = F)

################################
rm(list=ls())
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(scater)
library(SingleCellExperiment)
library(dplyr) 
library(ggplot2)
library(gridExtra)
library(cowplot)
setwd("/data5/zhangq/scRNA/MultiSample/NCTMerge9")
load("Kidney_NCT_merge.rda")
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/N_Anno')
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")
selectobject=ifnb.list[["Normal"]]

plot1<-DimPlot(selectobject, reduction = "umap", label = TRUE, pt.size = 0.4) + NoLegend()
plot2<-DimPlot(selectobject, reduction = "tsne", label = TRUE, pt.size = 0.4) + NoLegend()

p4=plot_grid(plot1,plot2,ncol = 2)
ggsave(p4,filename = "Kidney.integrated_Normal_Clusternum_umap_tsne.png",height = 6,width = 14)

head(selectobject@meta.data)
plott1=DimPlot(selectobject, reduction = "umap", label = TRUE, pt.size = 0.4,split.by = "orig.ident")
ggsave(plott1,filename = "Kidney.integrated_Normal_splitClusternum_umap.png",height = 6,width = 14)
plott2=DimPlot(selectobject, reduction = "tsne", label = TRUE, pt.size = 0.4,split.by = "orig.ident")
ggsave(plott2,filename = "Kidney.integrated_Normal_splitClusternum_tsne.png",height = 6,width =14)

cluster_cellnum=as.matrix(table(selectobject@meta.data$seurat_clusters))
cluster_cellnum=cbind(rownames(cluster_cellnum),cluster_cellnum)
colnames(cluster_cellnum)=c("Cluster number","Cell number")
write.table(cluster_cellnum,"Normal_cluster_cellnum.txt",col.names = T,row.names = F,sep = "\t",quote = F)

