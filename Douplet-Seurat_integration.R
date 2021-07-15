rm(list = ls())
options(future.globals.maxSize=160*1024^3)
library(Seurat)
library(Matrix)
library(dplyr)
library(DoubletFinder)
library(ggplot2)
options(stringsAsFactors=FALSE)
setwd("/your_path/")
dir.create("filter")
dir.create("marker")
samplesset = c("n0804","n0903","n1024","n1117"
              ,"c0804","c0826","c0903","c1020","c1024","c1102","c1117","c1202","c1205"
			  ,"t0804","t0826","t0903","t1020","t1024","t1102","t1117","t1202","t1205")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

result.matrix = matrix(0,length(samplesset)[1],5)
i = 0
Kidney.list = list()
for(sample_s in samplesset){
  i = i + 1
  print(sample_s)
  datadir = paste("/data5/zhangq/scRNA/",sample_s,"/outs/filtered_feature_bc_matrix/",sep="")
  Kidney.data.i <- Read10X(data.dir = datadir)
  result.matrix[i,1] = sample_s
  result.matrix[i,2] = dim(Kidney.data.i)[2]
  Kidney.seurat <- CreateSeuratObject(counts = Kidney.data.i, min.cells = 3, min.features = 200, project = sample_s)
  {
    if(substring(sample_s,1,1)=="n")
    {
      Kidney.seurat$group="Normal"
    }
    else if(substring(sample_s,1,1)=="c")
    {
      Kidney.seurat$group="Cancer"
    }
    else if(substring(sample_s,1,1)=="t")
    {
      Kidney.seurat$group="Thrombus"
    }
    
  }
  Kidney.seurat[["percent.mito"]] <- PercentageFeatureSet(Kidney.seurat, pattern = "^MT-")
  print(dim(Kidney.seurat))
  HB.genes_total=c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m=match(HB.genes_total,rownames(Kidney.seurat@assays$RNA))
  HB.genes=rownames(Kidney.seurat@assays$RNA)[HB_m]
  HB.genes=HB.genes[!is.na(HB.genes)]
  Kidney.seurat[["percent.HB"]]=PercentageFeatureSet(Kidney.seurat,features=HB.genes)
  head(Kidney.seurat@meta.data)[,c(2,3,4,5)]
  
  dpi = 300
  png(file=paste('filter/',sample_s,"_qc.png",sep=''), width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
  print(VlnPlot(object = Kidney.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()
  
  png(file=paste('filter/',sample_s,"_umi-mito.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  print(FeatureScatter(object = Kidney.seurat, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  dev.off()
  
  png(file=paste('filter/',sample_s,"_umi-gene.png",sep=''), width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
  print(FeatureScatter(object = Kidney.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  dev.off()
  
  result.matrix[i,3] = dim(Kidney.seurat)[2]
  
  Kidney.seurat <- subset(Kidney.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 50)
  dim(Kidney.seurat)
  Kidney.seurat <- NormalizeData(Kidney.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  Kidney.seurat <- CellCycleScoring(Kidney.seurat, g2m.features=g2m.genes, s.features=s.genes)
  Kidney.seurat <- FindVariableFeatures(Kidney.seurat, selection.method = "vst", nfeatures = 2000)
  result.matrix[i,4] = dim(Kidney.seurat)[2]
  Kidney.seurat=ScaleData(Kidney.seurat,verbose=FALSE)
  Kidney.seurat <- RunPCA(Kidney.seurat,verbose = F,features = VariableFeatures(Kidney.seurat)
                            ,npcs = 50)
  Kidney.seurat <- FindNeighbors(Kidney.seurat, dims = 1:50)
  Kidney.seurat <- FindClusters(Kidney.seurat)
  #Doublet Finder
  ## pK Identification (no ground-truth) 
  set.seed(123)
  sweep.res.list_kidney <- paramSweep_v3(Kidney.seurat, PCs = 1:50, sct = F)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  mpK<-as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
  # 找最佳 nExp
  ## ex: annotations <- Kidney.seurat@meta.data$ClusteringResults
  
  annotations <- Kidney.seurat@meta.data$seurat_clusters
  set.seed(123)
  homotypic.prop <- modelHomotypic(annotations)
  doublet_rate=(length(Cells(Kidney.seurat))*0.031)/4000

  nExp_poi <- round(doublet_rate*ncol(Kidney.seurat@assays$RNA@data)) ## Assuming doublet_rate% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # Find Doublet
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  set.seed(123)
  Kidney.seurat <- doubletFinder_v3(Kidney.seurat, PCs = 1:50, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  #pANN_0.25_mpK_nExp_poi
  set.seed(123)
  Kidney.seurat <- doubletFinder_v3(Kidney.seurat, PCs = 1:50, pN = 0.25, pK = mpK,nExp = nExp_poi.adj, sct = FALSE)
  
  Kidney.seurat@meta.data[,"DF_hi.lo"] <- Kidney.seurat@meta.data [, paste0("DF.classifications_0.25_",mpK, "_", nExp_poi.adj)]  
  Kidney.seurat@meta.data$DF_hi.lo[which(Kidney.seurat@meta.data$DF_hi.lo == "Doublet" & Kidney.seurat@meta.data[, paste0("DF.classifications_0.25_", mpK, "_", nExp_poi)] == "Singlet")] <- "Doublet_lo"
  Kidney.seurat@meta.data$DF_hi.lo[which(Kidney.seurat@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"

  Kidney.seurat=subset(Kidney.seurat,subset=DF_hi.lo=="Singlet")
  
  result.matrix[i,5] = dim(Kidney.seurat)[2]
  Kidney.list[sample_s]=Kidney.seurat
}

result.dataframe = as.data.frame(result.matrix)
colnames(result.dataframe) = c('sample','before','filter1',"filter2","drop_Doublet")
write.table(result.dataframe,file='filter/statistics_filter.txt',row.names = FALSE,quote = FALSE,sep='\t')

Kidney <- FindIntegrationAnchors(object.list = Kidney.list, dims = 1:50)
Kidney.integrated <- IntegrateData(anchorset = Kidney, dims = 1:50,features.to.integrate = rownames(Kidney))

###first generate data and scale data in RNA assay
DefaultAssay(Kidney.integrated) <- "RNA"
Kidney.integrated[['percent.mito']] <- PercentageFeatureSet(Kidney.integrated, pattern = "^MT-")
Kidney.integrated <- NormalizeData(object = Kidney.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
Kidney.integrated <- CellCycleScoring(Kidney.integrated, g2m.features=g2m.genes, s.features=s.genes)
Kidney.integrated <- FindVariableFeatures(object = Kidney.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)

#Kidney.integrated <- ScaleData(Kidney.integrated,features =all.genes,verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
all.genes <- rownames(Kidney.integrated)
length(rownames(Kidney.integrated)) #24132
Kidney.integrated<- ScaleData(Kidney.integrated, features = all.genes)
Kidney.integrated<- ScaleData(Kidney.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito","S.Score","G2M.Score"))

##change to integrated assay
DefaultAssay(Kidney.integrated) <- "integrated"
dpi = 300
png(file="qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
VlnPlot(object = Kidney.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

png(file="umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
FeatureScatter(object = Kidney.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Run the standard workflow for visualization and clustering
# Kidney.integrated <- ScaleData(Kidney.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
# ,features = VariableFeatures(Kidney.integrated)
Kidney.integrated<- ScaleData(Kidney.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA","percent.mito","S.Score","G2M.Score"))

set.seed(123)
Kidney.integrated <- RunPCA(Kidney.integrated,verbose = F
                            ,features = VariableFeatures(Kidney.integrated)
                            ,npcs = 50)

Kidney.integrated <- ProjectDim(object = Kidney.integrated)
png(file="pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = Kidney.integrated,ndims = 50)
dev.off()

save(Kidney.integrated,file = "Kidney.integrated.nonCluster.rda")
###cluster
set.seed(123)
Kidney.integrated <- FindNeighbors(object = Kidney.integrated, dims = 1:50)
Kidney.integrated <- FindClusters(object = Kidney.integrated, resolution = 0.8) 

###tsne and umap
set.seed(123)
Kidney.integrated <- RunTSNE(object = Kidney.integrated, dims = 1:50)
set.seed(123)
Kidney.integrated <- RunUMAP(Kidney.integrated, reduction = "pca", dims = 1:50)

save(Kidney.integrated,file = "Kidney.integrated.nonplot.rda")
png(file="tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'tsne',label = TRUE,pt.size = 0.4)
dev.off()
png(file="umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'umap',label = TRUE,pt.size = 0.4)
dev.off()

png(file="Phase_tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated,reduction = "tsne",pt.size = 0.4,group.by  ='Phase',cols=c("yellow","#458B74","#483D8B"))
dev.off()

png(file="Phase_umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated,reduction = 'umap',pt.size = 0.4,group.by  ='Phase',cols=c("yellow","#458B74","#483D8B"))
dev.off()


png(file="group.orig.ident_umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'umap',group.by = "orig.ident",pt.size = 0.4)
dev.off()
png(file="group.orig.ident_tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'tsne',group.by = "orig.ident",pt.size = 0.4)
dev.off()

png(file="split.orig.ident_umap.png", width = dpi*36, height = dpi*18, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'umap',label = TRUE, split.by = "orig.ident",pt.size = 0.4,ncol=6)
dev.off()
png(file="split.orig.ident_tsne.png", width = dpi*36, height = dpi*18, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'tsne',label = TRUE, split.by = "orig.ident",pt.size = 0.4,ncol=6)
dev.off()

png(file="group_umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'umap',label = FALSE, group.by = 'group',pt.size = 0.4)
dev.off()
png(file="group_tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'tsne',label = FALSE, group.by = 'group',pt.size = 0.4)
dev.off()

head(Kidney.integrated@meta.data)

nonname1=colnames(Kidney.integrated@meta.data)[grep("pANN",colnames(Kidney.integrated@meta.data))] #36
nonname2=colnames(Kidney.integrated@meta.data)[grep("DF.classifications",colnames(Kidney.integrated@meta.data))]

nonname=c(nonname1,nonname2)
Kidney.integrated@meta.data=Kidney.integrated@meta.data[,-which(names(Kidney.integrated@meta.data)%in%nonname)]

nol1=ncol(Kidney.integrated@meta.data)

Kidney.integrated@meta.data[,nol1+1]=substring(Kidney.integrated@meta.data[,1],2,5)
colnames(Kidney.integrated@meta.data)[nol1+1]="PatientID"

nol2=ncol(Kidney.integrated@meta.data)
Kidney.integrated@meta.data[which(Kidney.integrated@meta.data[,4]=="Normal"),nol2+1]="N"

Kidney.integrated@meta.data[which(Kidney.integrated@meta.data[,4]=="Cancer"|Kidney.integrated@meta.data[,4]=="Thrombus"),nol2+1]="CT"

colnames(Kidney.integrated@meta.data)[nol2+1]="NorCorT"


dpi=300
png(file="PatientID_umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'umap',label = FALSE, group.by = 'PatientID',pt.size = 0.4)
dev.off()
png(file="PatientID_tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = Kidney.integrated, reduction = 'tsne',label = FALSE, group.by = 'PatientID',pt.size = 0.4)
dev.off()


####zhao marker  yaoyong DefaultAssay(Kidney.integrated) <- "RNA"
DefaultAssay(Kidney.integrated) <- "RNA"
#DefaultAssay(Kidney.integrated) <- "integrated"

save(Kidney.integrated, file = "Kidney_NCT_merge.rda")

dpi = 300
png(file="feature.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
VlnPlot(object = Kidney.integrated, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()



