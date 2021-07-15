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
setwd("/data5/zhangq/scRNA/thrombus/myRef2")
load("Tumor_Merge_Reference.rda")
normal=Project.se
setwd("/data5/zhangq/scRNA/MultiSample/NCTMerge9")
load("Kidney_NCT_merge.rda")
dir.create('/data5/zhangq/scRNA/MultiSample/NCTMerge9/C_Anno')
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/C_Anno')
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")

selectobject=ifnb.list[["Cancer"]]
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
# Beta cell-related markers
library(scater)
library(SingleCellExperiment)
Kidney.integrated.sce=as(Kidney.integrated.se,"SingleCellExperiment")

save(pred.combined,Kidney.integrated.sce,Kidney.integrated.se,file = "Kidney.integrated_Cancer_pred.rda")


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
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/C_Anno')
load("Kidney.integrated_Cancer_pred.rda")
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")
selectobject=ifnb.list[["Cancer"]]

pred_mat=as.data.frame(cbind(rownames(pred.combined),pred.combined$labels))
colnames(pred_mat)=c("key","Annotation")

new.cluster.ids <- pred_mat[,2]
names(new.cluster.ids) <- levels(selectobject)

selectobject <- RenameIdents(selectobject, new.cluster.ids)


plot1<-DimPlot(selectobject, reduction = "umap", label = TRUE, pt.size = 0.4) + NoLegend()
plot2<-DimPlot(selectobject, reduction = "tsne", label = TRUE, pt.size = 0.4) + NoLegend()

p4=plot_grid(plot1,plot2,ncol = 2, labels = c("Science ref","Science ref"))
ggsave(p4,filename = "Kidney.integrated_Cancer_Anno_cellname_umap_tsne.png",height = 6,width = 14)

write.csv(pred_mat,"2018_2019Mergeref_Cancer_Anno.csv",row.names = F,quote = F)

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
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/C_Anno')
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")
selectobject=ifnb.list[["Cancer"]]
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


save(pred.combined,Kidney.integrated.sce,Kidney.integrated.se,file = "Kidney.integrated_BP_Cancer_pred.rda")

#######################################
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
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/C_Anno')
load("Kidney.integrated_BP_Cancer_pred.rda")
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")
selectobject=ifnb.list[["Cancer"]]

pred_mat=as.data.frame(cbind(rownames(pred.combined),pred.combined$labels))
colnames(pred_mat)=c("key","Annotation")

new.cluster.ids <- pred_mat[,2]
names(new.cluster.ids) <- levels(selectobject)

selectobject <- RenameIdents(selectobject, new.cluster.ids)
plot1<-DimPlot(selectobject, reduction = "umap", label = TRUE, pt.size = 0.4) + NoLegend()
plot2<-DimPlot(selectobject, reduction = "tsne", label = TRUE, pt.size = 0.4) + NoLegend()

p4=plot_grid(plot1,plot2,ncol = 2,labels = c("SingleR ref","SingleR ref"))
ggsave(p4,filename = "Kidney.integrated_BP_Cancer_Anno_cellname_umap_tsne.png",height = 6,width = 14)

write.csv(pred_mat,"BlueprintEncodeData_Cancer_Anno.csv",row.names = F,quote = F)

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
setwd('/data5/zhangq/scRNA/MultiSample/NCTMerge9/C_Anno')
ifnb.list <- SplitObject(Kidney.integrated, split.by ="group")
selectobject=ifnb.list[["Cancer"]]

plot1<-DimPlot(selectobject, reduction = "umap", label = TRUE, pt.size = 0.4) + NoLegend()
plot2<-DimPlot(selectobject, reduction = "tsne", label = TRUE, pt.size = 0.4) + NoLegend()

p4=plot_grid(plot1,plot2,ncol = 2)
ggsave(p4,filename = "Kidney.integrated_Cancer_Clusternum_umap_tsne.png",height = 6,width = 14)

