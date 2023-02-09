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

pbmc.marker<-FindAllMarkers(pbmc,only.pos = TRUE)
write.table(pbmc.marker,"C:/Users/xjmik/Desktop/GSE98638_Totalmarker.txt",sep = "\t")
remove(pbmc.marker)
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
for (i in levels(Idents(pbmc))) {
  pbmcs<-subset(pbmc,idents = i)
  Idents(pbmcs)<-pbmcs@meta.data$Bulksample
  pbmcs.markers<-FindMarkers(pbmcs,only.pos = TRUE,ident.1 = "T",ident.2 = "N")
  write.table(pbmcs.markers,file = paste0("c:/Users/xjmik/Desktop/GSE98638_T_N",i,".txt"),sep = "\t")
}
remove(i,pbmcs,pbmcs.markers)
for (i in levels(Idents(pbmc))[c(1:9,11)]) {
  pbmcs<-subset(pbmc,idents = i)
  Idents(pbmcs)<-pbmcs@meta.data$Bulksample
  pbmcs.markers<-FindMarkers(pbmcs,only.pos = TRUE,ident.1 = "T",ident.2 = "P")
  write.table(pbmcs.markers,file = paste0("c:/Users/xjmik/Desktop/GSE98638_T_P",i,".txt"),sep = "\t")
  remove(pbmcs.markers)
}
remove(i,pbmcs)
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
Idents(pbmc)<-pbmc@meta.data$Bulksample2
pbmc.marker2<-FindAllMarkers(pbmc,only.pos = TRUE)
write.table(pbmc.marker2,"C:/Users/xjmik/Desktop/GSE98638_BulkTotalmarkers.txt",sep = "\t")
remove(pbmc.marker2)
for (i in levels(Idents(pbmc))) {
  pbmcs<-subset(pbmc,idents = i)
  Idents(pbmcs)<-pbmcs@meta.data$Bulksample
  pbmcs.markers<-FindMarkers(pbmcs,only.pos = TRUE,ident.1 = "T",ident.2 = "N")
  write.table(pbmcs.markers,file = paste0("c:/Users/xjmik/Desktop/GSE98638_Bulk_T_N",i,".txt"),sep = "\t")
}
remove(i,pbmcs,pbmcs.markers)
for (i in levels(Idents(pbmc))) {
  pbmcs<-subset(pbmc,idents = i)
  Idents(pbmcs)<-pbmcs@meta.data$Bulksample
  pbmcs.markers<-FindMarkers(pbmcs,only.pos = TRUE,ident.1 = "T",ident.2 = "P")
  write.table(pbmcs.markers,file = paste0("c:/Users/xjmik/Desktop/GSE98638_Bulk_T_P",i,".txt"),sep = "\t")
  remove(pbmcs.markers)
}
remove(i,pbmcs)
pbmc<-RunLDA(pbmc,labels = pbmc$majorCluster)
pbmc<-RunLDA(pbmc,labels = pbmc$Bulksample2,reduction.name = "lda_bulksample2")
pbmc<-RunTSNE(pbmc,dims = 1:20)
pbmc<-RunICA(pbmc)
pbmc<-RunSLSI(pbmc,graph = "RNA_snn")
pbmc<-RunSPCA(pbmc,graph = "RNA_snn")
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "lda_umap",dims = 1:11)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "lda_tsne",dims = 1:11)
Idents(pbmc)<-pbmc@meta.data$majorCluster
DimPlot(pbmc,reduction = "lda_umap")
DimPlot(pbmc,reduction = "lda_tsne")
pbmc<-RunUMAP(pbmc,reduction = "lda_bulksample2",reduction.name = "lda_bulksample2_umap",dims = 1:2)
pbmc<-RunTSNE(pbmc,reduction = "lda_bulksample2",reduction.name = "lda_bulksample2_tsne",dims = 1:2)
Idents(pbmc)<-pbmc@meta.data$Bulksample2
DimPlot(pbmc,reduction = "lda_bulksample2_umap")
DimPlot(pbmc,reduction = "lda_bulksample2_tsne")
saveRDS(pbmc,file = "c:/Users/xjmik/Desktop/GSE98638/pbmc.rds")
DimPlot(pbmc,reduction = "lda_tsne")
Idents(pbmc)<-pbmc@meta.data$majorCluster
pbmc_avg<-AverageExpression(pbmc)
write.table(pbmc_avg$RNA,file = "c:/Users/xjmik/Downloads/GSE98638/average_Expression_cluster.txt",sep = "\t")
remove(pbmc_avg)
pbmc.markers<-FindMarkers(pbmc,only.pos = TRUE,ident.1 = "C08_CD4-CTLA4",ident.2 = "C07_CD4-FOXP3")
write.table(pbmc.markers,file = "c:/Users/xjmik/Desktop/GSE98638/Total_CTLA4_vs_FOXP3.txt",sep = "\t")
remove(pbmc.markers)
Idents(pbmc)<-pbmc@meta.data$Bulksample
DimPlot(pbmc,reduction = "lda_tsne")
FeaturePlot(pbmc,features = "CTLA4",reduction = "lda_tsne")
Idents(pbmc)<-pbmc@meta.data$majorCluster
DotPlot(pbmc, features = c("S100A4","DYNLL1","GPX1","ATP5E","MYL6")) + RotatedAxis()
venn<-read.table("c:/Users/xjmik/Desktop/GSE98638/Bulk3_venn_dotplot.txt",sep = "\t")
venn<-venn[,1]
DotPlot(pbmc, features = venn) + RotatedAxis()
remove(venn)
Idents(pbmc)<-pbmc$Bulksample
pbmcn<-subset(pbmc,idents = c("T","N","P"))
DotPlot(pbmcn, features = c("S100A4","DYNLL1","GPX1","ATP5E","MYL6")) + RotatedAxis()
venn<-read.table("c:/Users/xjmik/Desktop/GSE98638/Bulk3_venn_dotplot.txt",sep = "\t")
venn<-venn[,1]
DotPlot(pbmcn, features = venn) + RotatedAxis()
remove(venn,pbmcn)
Idents(pbmc)<-pbmc$majorCluster

h<-read.table("c:/Users/xjmik/Downloads/hallmarks.txt",sep = "\t")
rownames(h)<-h$V1
h<-as.data.frame(t(h))
h<-h[-1,]
h<-h[-1,]
h<-as.list(h)
library(GSVA)
gsva<-gsva(expr = as.matrix(pbmc@assays$RNA@data),gset.idx.list = h,kcdf="Poisson",parallel.sz=20,method ="ssgsea")
write.table(gsva,file="c:/Users/xjmik/Desktop/GSE98638/HallmarksGSVA.txt",sep="\t")

h<-read.table("c:/Users/xjmik/Downloads/c2CGP.txt",sep = "\t")
rownames(h)<-h$V1
h<-as.data.frame(t(h))
h<-h[-1,]
h<-h[-1,]
h<-as.list(h)
library(GSVA)
gsva_2<-gsva(expr = as.matrix(pbmc@assays$RNA@data),gset.idx.list = h,kcdf="Poisson",parallel.sz=20,method ="ssgsea")
write.table(gsva,file="c:/Users/xjmik/Desktop/GSE98638/c2CGPGSVA.txt",sep="\t")

h<-read.table("c:/Users/xjmik/Downloads/c2KEGG.txt",sep = "\t")
rownames(h)<-h$V1
h<-as.data.frame(t(h))
h<-h[-1,]
h<-h[-1,]
h<-as.list(h)
library(GSVA)
gsva_3<-gsva(expr = as.matrix(pbmc@assays$RNA@data),gset.idx.list = h,kcdf="Poisson",parallel.sz=20,method ="ssgsea")
write.table(gsva,file="c:/Users/xjmik/Desktop/GSE98638/c2KEGGGSVA.txt",sep="\t")

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

JMJD1Chi.data <- as.data.frame(Tregtumor1chi@assays$RNA@counts)
JMJD1Clo.data <- as.data.frame(Tregtumor1clo@assays$RNA@counts)
Tregpbmc.data <- as.data.frame(Tregpbmc@assays$RNA@counts)
remove(Tregtumor1chi,Tregtumor1clo,Tregpbmc,Tregtumor1cme)

for (i in 1:length(colnames(JMJD1Chi.data))) {
  colnames(JMJD1Chi.data)[i] <- paste(colnames(JMJD1Chi.data)[i],"JMJD1Chi",i,sep = "-")  
}

for (i in 1:length(colnames(JMJD1Clo.data))) {
  colnames(JMJD1Clo.data)[i] <- paste(colnames(JMJD1Clo.data)[i],"JMJD1Clo",i,sep = "-")  
}

for (i in 1:length(colnames(Tregpbmc.data))) {
  colnames(Tregpbmc.data)[i] <- paste(colnames(Tregpbmc.data)[i],"Tregpbmc",i,sep = "-")  
}

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
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",split.by = "group")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",split.by = "group")
pbmc <- RunLDA(pbmc,labels = pbmc@meta.data$seurat_clusters)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "lda_umap",dims = 1:4)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "lda_tsne",dims = 1:4)

# Idents(pbmc)<-pbmc@meta.data$group
# pbmc<-RunLDA(pbmc,labels = pbmc$group)
# DimPlot(pbmc, reduction = "lda")
# pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "lda_umap",dims = 1:2)
# pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "lda_tsne",dims = 1:2)

DimPlot(pbmc,reduction = "lda_umap")
DimPlot(pbmc,reduction = "lda_tsne")

# pbmc_new_3<-as.SingleCellExperiment(pbmc)
# 
# library(slingshot)
# sce_6<-slingshot(pbmc_new_3,clusterLabels = "group",reducedDim = 'LDA',start.clus = "pbmc",dist.method= "slingshot") 
# library(TrajectoryUtils)
# library(grDevices)
# library(RColorBrewer)
# library(ggplot2)
# library(ggsn)
# library(scales)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sce_6$slingPseudotime_1, breaks=100)]
# plot(reducedDims(sce_6)$LDA, col = plotcol, pch=16, asp = 1) 
# lines(SlingshotDataSet(sce_6), lwd=2, col='black')
# data<-sce_6@int_colData$reducedDims@listData$LDA
# color<-data.frame(sce_6$slingPseudotime_1,plotcol)
# rownames(color)<-rownames(data)
# data<-cbind(data,color)
# colnames(data)[3]<-"Pseudotime"
# colnames(data)[4]<-"color"
# remove(color)
# library(ggplot2)
# library(scales)
# library(ggsn)
# ggplot(data, aes(x = LDA_1, y = LDA_2,fill=Pseudotime)) +
#   geom_point(col = plotcol ) + 
#   theme_classic() +
#   scale_fill_gradientn(colors = colors)
# remove(plotcol,colors,all.genes,sce_6,pbmc_new_3,data)
remove(all.genes)
FeaturePlot(pbmc,features = "JMJD1C",reduction = "tsne")
FeaturePlot(pbmc,features = "JMJD1C",reduction = "tsne")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "JMJD1C",reduction = "tsne",min.cutoff = "q10",max.cutoff = "q90")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
DotPlot(pbmc, features = "JMJD1C") + RotatedAxis() + scale_size(range=c(5, 20))
VlnPlot(pbmc,features = "JMJD1C",sort = TRUE)
new.cluster.ids <- c("cluster4&3","cluster1&2&5&6","cluster0")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", pt.size = 0.5) 
DotPlot(pbmc, features = "JMJD1C") + RotatedAxis() + scale_size(range=c(5, 20))
VlnPlot(pbmc,features = "JMJD1C",sort = TRUE)
pbmc@meta.data$GROUP<-Idents(pbmc)
library(monocle)
trace('project2MST',edit = T,where = asNamespace("monocle"))
data<-as.sparse(pbmc@assays$RNA@counts)
pd <-pbmc@meta.data
pd <- new('AnnotatedDataFrame', data = pd)  
fData<-data.frame(gene_short_name = row.names(data),row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data, phenoData = pd,featureData = fd,expressionFamily=negbinomial())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
# 
# 
# 
pbmc.marker<-FindAllMarkers(pbmc,only.pos = TRUE,logfc.threshold = 0,min.pct = 0)
# diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = "~group",cores = 20)
# write.table(diff_test_res,"c:/Users/xjmik/Desktop/difftest.txt",sep = "\t")
write.table(pbmc.marker,"c:/Users/xjmik/Desktop/pbmcmarkers_1.txt",sep = "\t")
pbmcmarkers<-read.table("c:/Users/xjmik/Desktop/pbmcmarkers_1.txt",sep = "\t",header = TRUE,row.names = 1)
# diff_test_res<-read.table("c:/Users/xjmik/Desktop/difftest.txt",sep = "\t",header = TRUE,row.names = 1)
# remove(diff_test_res)
pbmcmarkers_new<-pbmcmarkers[which(pbmcmarkers$p_val_adj < 0.05 & pbmcmarkers$avg_log2FC > 0.25 & pbmcmarkers$pct.1 > 0.1 & pbmcmarkers$pct.2 > 0.1),]
pbmcmarkers_new$filter<-pbmcmarkers_new$pct.1 - pbmcmarkers_new$pct.2
pbmcmarkers_new_new<-pbmcmarkers_new[which(pbmcmarkers_new$filter > median(pbmcmarkers_new$filter)),]
ordering_genes <- pbmcmarkers_new_new$gene
monocle_cds <-setOrderingFilter(monocle_cds,ordering_genes = ordering_genes)
monocle_cds <-reduceDimension(monocle_cds,reduction_method = "DDRTree",max_components = 2)

monocle_cds <-orderCells(monocle_cds )
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "GROUP",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters",cell_size = 0.75)

remove(data,fd,fData,pbmcmarkers,pbmcmarkers_new,pbmcmarkers_new_new,pd,ordering_genes,new.cluster.ids)
remove(pbmc.marker)
# remove(data,diff_test_res,fd,fData,pbmc.marker,pd,new.cluster.ids)
#remove(pd,monocle_cds,fData,fd,data)
# blast_genes <- row.names(subset(fData(monocle_cds),
#                                 gene_short_name %in% c("JMJD1C")))
# plot_genes_jitter(monocle_cds[blast_genes,],
#                   grouping = "State",
#                   min_expr = 0.1)

my_genes <- row.names(subset(fData(monocle_cds),
                             gene_short_name %in% c("JMJD1C")))
cds_subset <- monocle_cds[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")
plotdf=pData(monocle_cds)
library(ggridges)
mycolor<-c("#619CFF","#00BA38","#F8766D")
ggplot(plotdf, aes(x=Pseudotime,y=GROUP,fill=GROUP))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+scale_fill_manual(values = mycolor)
# ggplot(plotdf, aes(x=Pseudotime,fill=GROUP,y=))+
#   geom_density_ridges(scale=1) +
#   geom_vline(xintercept = c(5,10),linetype=2)+
#   scale_y_discrete("")+
#   theme_minimal()+
#   theme(
#     panel.grid = element_blank()
#   )

remove(cds_subset,data,fd,fData,pbmcmarkers,pbmcmarkers_new,pbmcmarkers_new_new,pd,ordering_genes)
remove(blast_genes,my_genes)
remove(mycolor,plotdf)
Idents(pbmc)<-pbmc@meta.data$group

Tregtumor1chi<-subset(pbmc,idents = "hi")
Tregtumor1clo<-subset(pbmc,idents = "lo")
Tregpbmc<-subset(pbmc,idents = "pbmc")

JMJD1Chi.data <- as.data.frame(Tregtumor1chi@assays$RNA@counts)
JMJD1Clo.data <- as.data.frame(Tregtumor1clo@assays$RNA@counts)
Tregpbmc.data <- as.data.frame(Tregpbmc@assays$RNA@counts)
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

JMJD1Chi.metadata<-data.frame(colnames(JMJD1Chi.data),rep("hi",length(colnames(JMJD1Chi.data))))
JMJD1Clo.metadata<-data.frame(colnames(JMJD1Clo.data),rep("lo",length(colnames(JMJD1Clo.data))))
Tregpbmc.metadata<-data.frame(colnames(Tregpbmc.data),rep("pbmc",length(colnames(Tregpbmc.data))))
colnames(JMJD1Chi.metadata)<-c("barcode","group")
colnames(JMJD1Clo.metadata)<-c("barcode","group")
colnames(Tregpbmc.metadata)<-c("barcode","group")

rownames(JMJD1Chi.metadata)<-JMJD1Chi.metadata[,1]
rownames(JMJD1Clo.metadata)<-JMJD1Clo.metadata[,1]
rownames(Tregpbmc.metadata)<-Tregpbmc.metadata[,1]


JMJD1Chi <- CreateSeuratObject(counts = JMJD1Chi.data, project = "IMMUNE_JMJD1Chi",meta.data = JMJD1Chi.metadata)
JMJD1Chi$type <- "JMJD1Chi"
JMJD1Chi <- NormalizeData(JMJD1Chi, verbose = FALSE)
JMJD1Chi <- FindVariableFeatures(JMJD1Chi, selection.method = "vst", nfeatures = 2000)


JMJD1Clo <- CreateSeuratObject(counts = JMJD1Clo.data, project = "IMMUNE_JMJD1Clo",meta.data = JMJD1Clo.metadata)
JMJD1Clo$type <- "JMJD1Clo"
JMJD1Clo <- NormalizeData(JMJD1Clo, verbose = FALSE)
JMJD1Clo <- FindVariableFeatures(JMJD1Clo, selection.method = "vst", nfeatures = 2000)

Tregpbmc <- CreateSeuratObject(counts = Tregpbmc.data, project = "IMMUNE_Tregpbmc",meta.data = Tregpbmc.metadata)
Tregpbmc$type <- "Tregpbmc"
Tregpbmc <- NormalizeData(Tregpbmc, verbose = FALSE)
Tregpbmc <- FindVariableFeatures(Tregpbmc, selection.method = "vst", nfeatures = 2000)
remove(JMJD1Chi.data,JMJD1Clo.data,JMJD1Chi.metadata,JMJD1Clo.metadata,Tregpbmc.data,Tregpbmc.metadata)
remove(i)

immune.anchors <- FindIntegrationAnchors(object.list = list(JMJD1Chi,JMJD1Clo,Tregpbmc), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- ScaleData(immune.combined, verbose = FALSE,assay = "RNA")
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.8)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "type")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
p1
p2
DimPlot(immune.combined, reduction = "umap", split.by = "type")

immune.combined <- RunTSNE(immune.combined, dims = 1:20)
DimPlot(immune.combined, reduction = "tsne",split.by = "type")
DimPlot(immune.combined, reduction = "tsne",group.by  = "type")
DimPlot(immune.combined, reduction = "tsne")

immune.combined<-RunLDA(immune.combined,labels = immune.combined@meta.data$group,assay = "RNA",features = rownames(immune.combined))
immune.combined<-RunLDA(immune.combined,labels = immune.combined@meta.data$group,assay = "integrated",features = rownames(immune.combined),reduction.name = "LDA_integrated")
immune.combined<-RunTSNE(immune.combined,reduction = "lda",reduction.name = "lda_tsne",dims = 1:2)
immune.combined<-RunTSNE(immune.combined,reduction = "LDA_integrated",reduction.name = "lda_tsne_integrated",dims = 1:2)
Idents(immune.combined)<-immune.combined$group
DimPlot(immune.combined,reduction = "lda_tsne")
DimPlot(immune.combined,reduction = "lda_tsne_integrated")
remove(Tregpbmc,p1,p2,JMJD1Chi,JMJD1Clo,immune.anchors)
Treghilo<-subset(pbmc,idents = c("hi","lo"))
a<-as.data.frame(Treghilo@assays$RNA@data)
b<-a[(rowMeans(a[,c(1:102)])>0),]
c<-b[(rowMeans(b[,c(103:582)])>0),]
write.table(c,"c:/Users/xjmik/Desktop/hilo.txt",sep = "\t")
remove(a,b,c)
remove(Treghilo)
