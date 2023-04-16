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

CD4<-subset(pbmc,idents = c("C06_CD4-CCR7","C11_CD4-GNLY","C09_CD4-GZMA","C08_CD4-CTLA4","C07_CD4-FOXP3","C10_CD4-CXCL13"))
VlnPlot(CD4,features = "JMJD1C",sort = TRUE)
new.cluster.ids <- c("pbmcTC", "TumorTC", "TumorTC", "TumorTr", "pbmcTr", "TumorTC")
names(new.cluster.ids) <- levels(CD4)
CD4 <- RenameIdents(CD4, new.cluster.ids)
VlnPlot(CD4,features = "JMJD1C",sort = TRUE)
Treg<-subset(CD4,idents = c("TumorTr","pbmcTr"))
TC<-subset(CD4,idents = c("TumorTC","pbmcTC"))
library(ggpubr)
VlnPlot(Treg,features = "JMJD1C",sort = TRUE)+stat_compare_means()
VlnPlot(TC,features = "JMJD1C",sort = TRUE)+stat_compare_means()
