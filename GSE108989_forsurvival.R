library(dplyr)
library(Seurat)

pbmc.data<-read.table("c:/Users/xjmik/Downloads/GSE108989/GSE108989_CRC.TCell.S11138.count.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
pbmc.data<-pbmc.data[!duplicated(pbmc.data$symbol),]
pbmc.data<-na.omit(pbmc.data)
rownames(pbmc.data)<-pbmc.data$symbol
pbmc.data<-pbmc.data[,-1]
pbmc.metadata<-read.table("c:/Users/xjmik/Downloads/GSE108989/metadata.txt",sep = "\t",header = TRUE,row.names = 1)
pbmc.metadata<-pbmc.metadata[colnames(pbmc.data),]
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
remove(all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

Idents(pbmc)<-pbmc@meta.data$majorCluster


for (i in 1:length(rownames(pbmc@meta.data))) {
  if(pbmc@meta.data$sampleType[i] == "PTH"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "PTC"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "PTR"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "TTH"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "TTR"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "TTC"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "NTY"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "TTY"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "PTY"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "NTH"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "NTR"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "NTC"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "NP7"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "PP7"){
    pbmc@meta.data$Bulksample[i] <- "P"
  }
  if(pbmc@meta.data$sampleType[i] == "TP7"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
}
remove(i)

for (i in 1:length(rownames(pbmc@meta.data))) {
  if(pbmc@meta.data$sampleType[i] == "PTH"){
    pbmc@meta.data$Bulksample2[i] <- "TH"
  }
  if(pbmc@meta.data$sampleType[i] == "PTC"){
    pbmc@meta.data$Bulksample2[i] <- "TC"
  }
  if(pbmc@meta.data$sampleType[i] == "PTR"){
    pbmc@meta.data$Bulksample2[i] <- "TR"
  }
  if(pbmc@meta.data$sampleType[i] == "TTH"){
    pbmc@meta.data$Bulksample2[i] <- "TH"
  }
  if(pbmc@meta.data$sampleType[i] == "TTR"){
    pbmc@meta.data$Bulksample2[i] <- "TR"
  }
  if(pbmc@meta.data$sampleType[i] == "TTC"){
    pbmc@meta.data$Bulksample2[i] <- "TC"
  }
  if(pbmc@meta.data$sampleType[i] == "NTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
  }
  if(pbmc@meta.data$sampleType[i] == "PTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
  }
  if(pbmc@meta.data$sampleType[i] == "TTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
  }
  if(pbmc@meta.data$sampleType[i] == "NTH"){
    pbmc@meta.data$Bulksample2[i] <- "TH"
  }
  if(pbmc@meta.data$sampleType[i] == "NTR"){
    pbmc@meta.data$Bulksample2[i] <- "TR"
  }
  if(pbmc@meta.data$sampleType[i] == "NTC"){
    pbmc@meta.data$Bulksample2[i] <- "TC"
  }
  if(pbmc@meta.data$sampleType[i] == "NP7"){
    pbmc@meta.data$Bulksample2[i] <- "P7"
  }
  if(pbmc@meta.data$sampleType[i] == "PP7"){
    pbmc@meta.data$Bulksample2[i] <- "P7"
  }
  if(pbmc@meta.data$sampleType[i] == "TP7"){
    pbmc@meta.data$Bulksample2[i] <- "P7"
  }
}
remove(i)


Idents(pbmc)<-pbmc@meta.data$majorCluster

Tregtumor<-subset(pbmc,idents = c("CD4_C12-CTLA4"))
Tregtumor <- NormalizeData(Tregtumor)
Tregtumor <- FindVariableFeatures(Tregtumor, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tregtumor)
Tregtumor <- ScaleData(Tregtumor, features = all.genes)
remove(all.genes)
Tregtumor <- RunPCA(Tregtumor, features = VariableFeatures(object = Tregtumor),)
Tregtumor <- RunUMAP(Tregtumor,dims = 1:30)
Tregtumor <- RunTSNE(Tregtumor,dims = 1:30)
Tregtumor <- FindNeighbors(Tregtumor, dims = 1:7,k.param = 10,n.trees = 10,reduction = "pca",annoy.metric = "euclidean")
Tregtumor <- FindClusters(Tregtumor, resolution = 0.3,algorithm = 2,modularity.fxn = 1)
Tregtumor.marker<-FindAllMarkers(Tregtumor,only.pos = TRUE)
DotPlot(Tregtumor, features = "JMJD1C") + RotatedAxis()
VlnPlot(Tregtumor,features = "JMJD1C",sort = TRUE)

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
Tregtumor<-RunLDA(Tregtumor,labels = Tregtumor$seurat_clusters)
Tregtumor<-RunUMAP(Tregtumor,reduction = "lda",reduction.name = "lda_umap",dims = 1:3)
Tregtumor<-RunTSNE(Tregtumor,reduction = "lda",reduction.name = "lda_tsne",dims = 1:3)
DimPlot(Tregtumor,reduction = "lda_umap")
DimPlot(Tregtumor,reduction = "lda_tsne")
FeaturePlot(Tregtumor,features = "JMJD1C",label = TRUE,reduction = "lda_umap")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(Tregtumor,features = "JMJD1C",label = TRUE,reduction = "lda_tsne")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
new.cluster.ids <- c("lo","lo","lo","hi","hi","lo")
names(new.cluster.ids) <- levels(Tregtumor)
Tregtumor <- RenameIdents(Tregtumor, new.cluster.ids)
DimPlot(Tregtumor, reduction = "lda_umap", label = TRUE, pt.size = 0.5) + NoLegend()

remove(new.cluster.ids,Tregtumor.marker)
library(ggpubr)
VlnPlot(Tregtumor,features = "JMJD1C",sort = TRUE)+stat_compare_means()
Tregtumor1chi<-subset(Tregtumor,idents = "hi")
Tregtumor1clo<-subset(Tregtumor,idents = "lo")

Tregpbmc<-subset(pbmc,idents = "CD4_C10-FOXP3")
remove(Tregtumor)
JMJD1Chi.data <- as.data.frame(Tregtumor1chi@assays$RNA@counts)
JMJD1Clo.data <- as.data.frame(Tregtumor1clo@assays$RNA@counts)
Tregpbmc.data <- as.data.frame(Tregpbmc@assays$RNA@counts)
remove(Tregtumor1chi,Tregtumor1clo,Tregpbmc)


JMJD1Chi.metadata<-data.frame(colnames(JMJD1Chi.data),rep("hi",length(colnames(JMJD1Chi.data))))
JMJD1Clo.metadata<-data.frame(colnames(JMJD1Clo.data),rep("lo",length(colnames(JMJD1Clo.data))))
Tregpbmc.metadata<-data.frame(colnames(Tregpbmc.data),rep("pbmc",length(colnames(Tregpbmc.data))))
colnames(JMJD1Chi.metadata)<-c("barcode","group")
colnames(JMJD1Clo.metadata)<-c("barcode","group")
colnames(Tregpbmc.metadata)<-c("barcode","group")
pbmc.metadata<-rbind(JMJD1Chi.metadata,JMJD1Clo.metadata,Tregpbmc.metadata)
rownames(pbmc.metadata)<-pbmc.metadata[,1]
pbmc.data<-cbind(JMJD1Chi.data,JMJD1Clo.data,Tregpbmc.data)
remove(i,JMJD1Chi.metadata,JMJD1Clo.metadata,Tregpbmc.metadata,JMJD1Chi.data,JMJD1Clo.data,Tregpbmc.data)
pbmc_raw<-pbmc
remove(pbmc)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",split.by = "group")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",split.by = "group")
pbmc <- RunLDA(pbmc,labels = pbmc@meta.data$seurat_clusters)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "lda_umap",dims = 1:9)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "lda_tsne",dims = 1:9)

pbmc<-RunLDA(pbmc,labels = pbmc$group,reduction.name = "lda2")
Idents(pbmc)<-pbmc$group
DimPlot(pbmc, reduction = "lda2")
pbmc<-RunUMAP(pbmc,reduction = "lda2",reduction.name = "lda2_umap",dims = 1:2)
pbmc<-RunTSNE(pbmc,reduction = "lda2",reduction.name = "lda2_tsne",dims = 1:2)

DimPlot(pbmc,reduction = "lda2_umap")
DimPlot(pbmc,reduction = "lda2_tsne")
DimPlot(pbmc, reduction = "lda2_umap",split.by = "group")
DimPlot(pbmc, reduction = "lda2_tsne",split.by = "group")
remove(all.genes)

hilomarkers<-FindMarkers(pbmc,only.pos = TRUE,ident.1 = "hi",ident.2 = "lo")
write.table(hilomarkers,"c:/Users/xjmik/Desktop/GSE108989hilomarkers.txt",sep = "\t")
remove(hilomarkers)
name<-colnames(subset(pbmc,idents = "hi"))

pbmc_raw$GROUP<-rep("other",length(colnames(pbmc_raw)))
for (i in 1:length(name)) {
  for (j in 1:length(colnames(pbmc_raw))) {
    if(name[i] == colnames(pbmc_raw)[j]){
      pbmc_raw@meta.data$GROUP[j] <- "hi"
    }
  }
}

remove(i,j)
Idents(pbmc_raw)<-pbmc_raw$GROUP
markers1<-FindMarkers(pbmc_raw,ident.1 = "hi",ident.2 = "other",only.pos = TRUE)
write.table(markers1,file = "c:/Users/xjmik/Desktop/GSE108989hivsother.txt",sep = "\t")
remove(name)
remove(markers1)
