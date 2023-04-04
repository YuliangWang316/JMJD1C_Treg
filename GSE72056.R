library(dplyr)
library(Seurat)
library(patchwork)
pbmc.data <- read.table("c:/Users/xjmik/Downloads/GSE72056_melanoma_single_cell_revised_v2.txt",sep = "\t",header = TRUE)
pbmc.data<-pbmc.data[-1,]
pbmc.data<-pbmc.data[-1,]
pbmc.data<-pbmc.data[-1,]
pbmc.data<-pbmc.data[!duplicated(pbmc.data$Cell),]
rownames(pbmc.data)<-pbmc.data$Cell
pbmc.data<-pbmc.data[,-1]
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:30)
DimPlot(pbmc, reduction = "umap")
a<-read.table("c:/Users/xjmik/Downloads/Treg_1C_signature.txt",sep = "\t",header = TRUE)
a_new<-a[1:length(rownames(a)),]
pbmc<-AddModuleScore(pbmc,features = a_new,name = "Treg_1C")
d<-pbmc@meta.data
Treg_1C<-dplyr::select(d,starts_with("Treg_1C"))
pbmc@meta.data$Treg_1C<-apply(Treg_1C,MARGIN = 1,median)
FeaturePlot(pbmc,features = "Treg_1C")
library(RColorBrewer)
library(ggpubr)
FeaturePlot(pbmc,features = "Treg_1C",reduction = "umap",label = FALSE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "FOXP3",reduction = "umap",label = FALSE)+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
