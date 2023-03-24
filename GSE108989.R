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
Tregtumor <- FindNeighbors(Tregtumor, dims = 1:5,k.param = 4,n.trees = 200,reduction = "pca",annoy.metric = "euclidean")
Tregtumor <- FindClusters(Tregtumor, resolution = 0.4,algorithm = 1,modularity.fxn = 1)
Tregtumor.marker<-FindAllMarkers(Tregtumor,only.pos = TRUE)
DotPlot(Tregtumor, features = "JMJD1C") + RotatedAxis()
VlnPlot(Tregtumor,features = "JMJD1C",sort = TRUE)
