library(dplyr)
library(Seurat)

pbmc.data<-read.table("c:/Users/xjmik/Downloads/GSE98638/GSE98638_HCC.TCell.S5063.count.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
pbmc.data<-pbmc.data[!duplicated(pbmc.data$symbol),]
pbmc.data<-na.omit(pbmc.data)
rownames(pbmc.data)<-pbmc.data$symbol
pbmc.data<-pbmc.data[,-1]
pbmc.metadata<-read.table("c:/Users/xjmik/Downloads/GSE98638/metadata.txt",sep = "\t",header = TRUE,row.names = 1)
pbmc.metadata<-pbmc.metadata[colnames(pbmc.data),]
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
remove(all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc,dims = 1:20)
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
  if(pbmc@meta.data$sampleType[i] == "JTH"){
    pbmc@meta.data$Bulksample[i] <- "J"
  }
  if(pbmc@meta.data$sampleType[i] == "JTR"){
    pbmc@meta.data$Bulksample[i] <- "J"
  }
  if(pbmc@meta.data$sampleType[i] == "JTC"){
    pbmc@meta.data$Bulksample[i] <- "J"
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
  if(pbmc@meta.data$sampleType[i] == "JTH"){
    pbmc@meta.data$Bulksample2[i] <- "TH"
  }
  if(pbmc@meta.data$sampleType[i] == "JTR"){
    pbmc@meta.data$Bulksample2[i] <- "TR"
  }
  if(pbmc@meta.data$sampleType[i] == "JTC"){
    pbmc@meta.data$Bulksample2[i] <- "TC"
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
}
remove(i)

Idents(pbmc)<-pbmc$majorCluster
Treg<-subset(pbmc,idents = c("C08_CD4-CTLA4"))
Treg <- NormalizeData(Treg)
Treg <- FindVariableFeatures(Treg, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Treg)
Treg <- ScaleData(Treg, features = all.genes)
Treg <- RunPCA(Treg, features = VariableFeatures(object = Treg))
ElbowPlot(Treg)
Treg <- FindNeighbors(Treg, dims = 1:10)
Treg <- FindClusters(Treg, resolution = 0.4)
Treg <- RunUMAP(Treg, dims = 1:10)
Treg <- RunTSNE(Treg, dims = 1:10)
VlnPlot(Treg,features = "JMJD1C",sort = TRUE)
FeaturePlot(Treg,features = "JMJD1C")
Tregmarkers<-FindAllMarkers(Treg,only.pos = TRUE)
DimPlot(Treg,reduction = "tsne")
DotPlot(Treg, features = "JMJD1C") + RotatedAxis()
Treg<-RunLDA(Treg,labels = Treg@meta.data$seurat_clusters)
DimPlot(Treg,reduction = "lda")

Treg<-RunUMAP(Treg,reduction = "lda",reduction.name = "lda_umap",dims = 1:4)
Treg<-RunTSNE(Treg,reduction = "lda",reduction.name = "lda_tsne",dims = 1:4)
DimPlot(Treg,reduction = "lda_umap")
DimPlot(Treg,reduction = "lda_tsne")

library(ggpubr)
library(RColorBrewer)
FeaturePlot(Treg,features = "JMJD1C",label = TRUE,reduction = "lda_tsne")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
remove(all.genes)
new.cluster.ids <- c("lo","lo","hi","hi","lo","lo")
names(new.cluster.ids) <- levels(Treg)
Treg <- RenameIdents(Treg, new.cluster.ids)
DimPlot(Treg, reduction = "lda_tsne", label = TRUE, pt.size = 0.5) + NoLegend()
VlnPlot(Treg,features = "JMJD1C",sort = TRUE)
remove(new.cluster.ids,Tregmarkers)

Tregtumor1chi<-subset(Treg,idents = "hi")
Tregtumor1clo<-subset(Treg,idents = "lo")
remove(Treg)
Tregpbmc<-subset(pbmc,idents = "C07_CD4-FOXP3")
CD8<-subset(pbmc,idents = c("C03_CD8-SLC4A10","C02_CD8-CX3CR1","C01_CD8-LEF1","C04_CD8-LAYN","C05_CD8-GZMK"))

JMJD1Chi.data <- as.data.frame(Tregtumor1chi@assays$RNA@counts)
JMJD1Clo.data <- as.data.frame(Tregtumor1clo@assays$RNA@counts)
Tregpbmc.data <- as.data.frame(Tregpbmc@assays$RNA@counts)
CD8.data<-as.data.frame(CD8@assays$RNA@counts)
remove(Tregtumor1chi,Tregtumor1clo,Tregpbmc)

for (i in 1:length(colnames(JMJD1Chi.data))) {
  colnames(JMJD1Chi.data)[i] <- paste(colnames(JMJD1Chi.data)[i],"JMJD1Chi",i,sep = "-")  
}

for (i in 1:length(colnames(JMJD1Clo.data))) {
  colnames(JMJD1Clo.data)[i] <- paste(colnames(JMJD1Clo.data)[i],"JMJD1Clo",i,sep = "-")  
}

for (i in 1:length(colnames(Tregpbmc.data))) {
  colnames(Tregpbmc.data)[i] <- paste(colnames(Tregpbmc.data)[i],"Tregpbmc",i,sep = "-")  
}

for (i in 1:length(colnames(CD8.data))) {
  colnames(CD8.data)[i] <- paste(colnames(CD8.data)[i],"CD8",i,sep = "-")  
}

JMJD1Chi.metadata<-data.frame(colnames(JMJD1Chi.data),rep("hi",length(colnames(JMJD1Chi.data))))
JMJD1Clo.metadata<-data.frame(colnames(JMJD1Clo.data),rep("lo",length(colnames(JMJD1Clo.data))))
Tregpbmc.metadata<-data.frame(colnames(Tregpbmc.data),rep("pbmc",length(colnames(Tregpbmc.data))))
CD8.metadata<-data.frame(colnames(CD8.data),CD8$majorCluster)

colnames(JMJD1Chi.metadata)<-c("barcode","group")
colnames(JMJD1Clo.metadata)<-c("barcode","group")
colnames(Tregpbmc.metadata)<-c("barcode","group")
colnames(CD8.metadata)<-c("barcode","group")
remove(CD8)

pbmc.metadata<-rbind(JMJD1Chi.metadata,JMJD1Clo.metadata,Tregpbmc.metadata,CD8.metadata)
rownames(pbmc.metadata)<-pbmc.metadata[,1]
pbmc.data<-cbind(JMJD1Chi.data,JMJD1Clo.data,Tregpbmc.data,CD8.data)
remove(i,JMJD1Chi.metadata,JMJD1Clo.metadata,Tregpbmc.metadata,JMJD1Chi.data,JMJD1Clo.data,Tregpbmc.data,CD8.data,CD8.metadata)
pbmc_raw<-pbmc
remove(pbmc)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",split.by = "group")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",split.by = "group")
pbmc <- RunLDA(pbmc,labels = pbmc@meta.data$group)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "lda_umap",dims = 1:7)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "lda_tsne",dims = 1:7)

Idents(pbmc)<-pbmc$group
DimPlot(pbmc,reduction = "lda_umap")
DimPlot(pbmc,reduction = "lda_tsne")


remove(all.genes)
DotPlot(pbmc, features = "JMJD1C") + RotatedAxis() + scale_size(range=c(5, 20))
VlnPlot(pbmc,features = "JMJD1C",sort = TRUE)

FeaturePlot(pbmc,features = "JMJD1C",reduction = "lda_umap")
FeaturePlot(pbmc,features = "JMJD1C",reduction = "lda_umap")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "JMJD1C",reduction = "lda_tsne")
FeaturePlot(pbmc,features = "JMJD1C",reduction = "lda_tsne")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))

