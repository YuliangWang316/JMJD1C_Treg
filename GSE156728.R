library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD4/byDataset/")
count<-data.frame()
metadata<-data.frame()
gene<-character()
c<-readRDS(a[1])
count<-c@assays$data@listData$norm_exprs
count<-as.matrix(count)
count<-as.data.frame(count)
metadata<-as.data.frame(c@colData@listData)
count<-count[,metadata$cellID]
colnames(count)<-metadata$cellID.uniq
gene<-rownames(count)
for (i in c(2:23,25:length(a))) {
  b<-readRDS(a[i])
  count_B<-b@assays$data@listData$norm_exprs
  count_B<-as.matrix(count_B)
  count_B<-as.data.frame(count_B)
  count_B_metadata<-as.data.frame(b@colData@listData)
  count_B<-count_B[,count_B_metadata$cellID]
  colnames(count_B)<-count_B_metadata$cellID.uniq
  genename<-rownames(count_B)
  d<-str_sub(genename[1],1,4)
  if(d == "ENSG"){
    e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
    e<-e[!duplicated(e$SYMBOL),]
    count_B<-count_B[e$ENSEMBL,]
    rownames(count_B)<-e$SYMBOL
    genename<-e$SYMBOL
  }
  f<-sort(genename)
  if(f[1] == "1"){
    g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
    g<-g[!duplicated(g$SYMBOL),]
    count_B<-count_B[g$ENTREZID,]
    rownames(count_B)<-g$SYMBOL
    genename<-g$SYMBOL
  }
  gene<-intersect(gene,genename)
  count_B<-count_B[gene,]
  count<-count[gene,]
  count<-cbind(count,count_B)
  metadata<-rbind(metadata,count_B_metadata)
  remove(count_B,count_B_metadata,b,d)
}
b<-readRDS(a[24])
count_B<-b@assays@data$norm_exprs
count_B<-as.matrix(count_B)
count_B<-as.data.frame(count_B)
count_B_metadata<-as.data.frame(b@colData@listData)
count_B<-count_B[,count_B_metadata$cellID]
colnames(count_B)<-count_B_metadata$cellID.uniq
genename<-rownames(count_B)
d<-str_sub(genename[1],1,4)
if(d == "ENSG"){
  e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
  e<-e[!duplicated(e$SYMBOL),]
  count_B<-count_B[e$ENSEMBL,]
  rownames(count_B)<-e$SYMBOL
  genename<-e$SYMBOL
}
f<-sort(genename)
if(f[1] == "1"){
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  count_B<-count_B[g$ENTREZID,]
  rownames(count_B)<-g$SYMBOL
  genename<-g$SYMBOL
}
gene<-intersect(gene,genename)
count_B<-count_B[gene,]
count<-count[gene,]
count<-cbind(count,count_B)
metadata<-rbind(metadata,count_B_metadata)
remove(count_B,count_B_metadata,b,d)
remove(c,e,g,a,f,gene,genename,i)
count_CD4<-count
metadata_CD4<-metadata
remove(count,metadata)

library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
setwd("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
a<-list.files("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/byDataset/")
count<-data.frame()
metadata<-data.frame()
gene<-character()
c<-readRDS(a[1])
count<-c@assays$data@listData$norm_exprs
count<-as.matrix(count)
count<-as.data.frame(count)
metadata<-as.data.frame(c@colData@listData)
count<-count[,metadata$cellID]
colnames(count)<-metadata$cellID.uniq
gene<-rownames(count)
for (i in c(2:23,25:length(a))) {
  b<-readRDS(a[i])
  count_B<-b@assays$data@listData$norm_exprs
  count_B<-as.matrix(count_B)
  count_B<-as.data.frame(count_B)
  count_B_metadata<-as.data.frame(b@colData@listData)
  count_B<-count_B[,count_B_metadata$cellID]
  colnames(count_B)<-count_B_metadata$cellID.uniq
  genename<-rownames(count_B)
  d<-str_sub(genename[1],1,4)
  if(d == "ENSG"){
    e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
    e<-e[!duplicated(e$SYMBOL),]
    count_B<-count_B[e$ENSEMBL,]
    rownames(count_B)<-e$SYMBOL
    genename<-e$SYMBOL
  }
  f<-sort(genename)
  if(f[1] == "1"){
    g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
    g<-g[!duplicated(g$SYMBOL),]
    count_B<-count_B[g$ENTREZID,]
    rownames(count_B)<-g$SYMBOL
    genename<-g$SYMBOL
  }
  gene<-intersect(gene,genename)
  count_B<-count_B[gene,]
  count<-count[gene,]
  count<-cbind(count,count_B)
  metadata<-rbind(metadata,count_B_metadata)
  remove(count_B,count_B_metadata,b,d)
}
b<-readRDS(a[24])
count_B<-b@assays@data$norm_exprs
count_B<-as.matrix(count_B)
count_B<-as.data.frame(count_B)
count_B_metadata<-as.data.frame(b@colData@listData)
count_B<-count_B[,count_B_metadata$cellID]
colnames(count_B)<-count_B_metadata$cellID.uniq
genename<-rownames(count_B)
d<-str_sub(genename[1],1,4)
if(d == "ENSG"){
  e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)  
  e<-e[!duplicated(e$SYMBOL),]
  count_B<-count_B[e$ENSEMBL,]
  rownames(count_B)<-e$SYMBOL
  genename<-e$SYMBOL
}
f<-sort(genename)
if(f[1] == "1"){
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  count_B<-count_B[g$ENTREZID,]
  rownames(count_B)<-g$SYMBOL
  genename<-g$SYMBOL
}
gene<-intersect(gene,genename)
count_B<-count_B[gene,]
count<-count[gene,]
count<-cbind(count,count_B)
metadata<-rbind(metadata,count_B_metadata)
remove(count_B,count_B_metadata,b,d)
remove(c,e,g,a,f,gene,genename,i)
count_CD8<-count
metadata_CD8<-metadata
remove(count,metadata)

metadata<-rbind(metadata_CD4,metadata_CD8)
count<-cbind(count_CD4,count_CD8)
e<-readRDS("C:/Users/xjmik/Downloads/data.expression/data/expression/CD4/integration/int.CD4.S35.meta.tb.rds")
f<-readRDS("c:/Users/xjmik/Downloads/data.expression/data/expression/CD8/integration/int.CD8.S35.meta.tb.rds")
g<-rbind(e,f)

remove(count_CD4,count_CD8,metadata_CD4,metadata_CD8,e,f)
h<-intersect(g$cellID.uniq,metadata$cellID.uniq)
g<-as.data.frame(g)
rownames(g)<-g$cellID.uniq
rownames(metadata)<-metadata$cellID.uniq
g_new<-g[h,]
metadata_new<-metadata[h,]
remove(metadata,g,h)
metadata<-cbind(g_new,metadata_new)
remove(g_new,metadata_new)
metadata<-metadata[colnames(count),]
library(Seurat)
pbmc<-CreateSeuratObject(counts = count,meta.data = metadata)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
remove(metadata,count)
#pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:30)
Idents(pbmc)<-pbmc@meta.data$meta.cluster

Treg<-subset(pbmc,ident = c("CD4.c18.Treg.RTKN2","CD4.c19.Treg.S1PR1","CD4.c20.Treg.TNFRSF9","CD4.c21.Treg.OAS1"))
Idents(Treg)<-Treg@meta.data$loc
Treg_T<-subset(Treg,ident = "T")
Treg_T <- NormalizeData(Treg_T)
Treg_T <- FindVariableFeatures(Treg_T, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Treg_T)

Treg_T <- ScaleData(Treg_T, features = all.genes)
Treg_T <- RunPCA(Treg_T, features = VariableFeatures(object = Treg_T))
ElbowPlot(Treg_T,ndims = 50)
Treg_T <- FindNeighbors(Treg_T, dims = 1:10)
Treg_T <- FindClusters(Treg_T, resolution = 0.4)
Treg_T <- RunUMAP(Treg_T, dims = 1:10)
Treg_Tmarkers<-FindAllMarkers(Treg_T,only.pos = TRUE)
JMJD1Chi<-subset(Treg_T,ident = c("1","5","7"))
DimPlot(Treg_T,label = TRUE)
Treg_T <- RunTSNE(Treg_T, dims = 1:10)
Treg_T <- RunLDA(Treg_T, labels = Treg_T$seurat_clusters)
remove(all.genes)
Treg_T<-RunUMAP(Treg_T,reduction = "lda",reduction.name = "lda_umap",dims = 1:9)
Treg_T<-RunTSNE(Treg_T,reduction = "lda",reduction.name = "lda_tsne",dims = 1:9)
DimPlot(Treg_T,label = TRUE,reduction = "tsne")
DimPlot(Treg_T,label = TRUE,reduction = "lda_umap")
DimPlot(Treg_T,label = TRUE,reduction = "lda_tsne")
VlnPlot(Treg_T,features = "JMJD1C",sort = TRUE,pt.size = 0)
library(ggplot2)
DotPlot(Treg_T, features = "JMJD1C") + RotatedAxis() + scale_size(range=c(5, 20))
library(ggpubr)
library(RColorBrewer)
FeaturePlot(Treg_T,features = "JMJD1C",label = TRUE,reduction = "umap")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(Treg_T,features = "JMJD1C",label = TRUE,reduction = "umap")
remove(JMJD1Chi)
new.cluster.ids <- c("lo","hi","lo","lo","lo","hi","lo","hi","lo","lo")
names(new.cluster.ids) <- levels(Treg_T)
Treg_T <- RenameIdents(Treg_T, new.cluster.ids)
DimPlot(Treg_T, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
VlnPlot(Treg_T,features = "JMJD1C",sort = TRUE,pt.size = 0)
remove(new.cluster.ids,Treg_Tmarkers)

Tregtumor1chi<-subset(Treg_T,idents = "hi")
Tregtumor1clo<-subset(Treg_T,idents = "lo")
remove(Treg_T)
Tregpbmc<-subset(Treg,idents = "P")
remove(Treg)
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
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.4)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",split.by = "group")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",split.by = "group")
pbmc <- RunLDA(pbmc,labels = pbmc@meta.data$seurat_clusters)
pbmc<-RunUMAP(pbmc,reduction = "lda",reduction.name = "lda_umap",dims = 1:10)
pbmc<-RunTSNE(pbmc,reduction = "lda",reduction.name = "lda_tsne",dims = 1:10)


DimPlot(pbmc,reduction = "lda_umap")
DimPlot(pbmc,reduction = "lda_tsne")



FeaturePlot(pbmc,features = "JMJD1C",reduction = "tsne")
FeaturePlot(pbmc,features = "JMJD1C",reduction = "tsne")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
FeaturePlot(pbmc,features = "JMJD1C",reduction = "tsne",min.cutoff = "q10",max.cutoff = "q90")+ scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(10))
DotPlot(pbmc, features = "JMJD1C") + RotatedAxis() + scale_size(range=c(5, 20))
VlnPlot(pbmc,features = "JMJD1C",sort = TRUE)
# new.cluster.ids <- c("cluster4&3","cluster1&2&5&6","cluster0")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "tsne", pt.size = 0.5) 
# DotPlot(pbmc, features = "JMJD1C") + RotatedAxis() + scale_size(range=c(5, 20))
# VlnPlot(pbmc,features = "JMJD1C",sort = TRUE)
remove(all.genes)
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
Idents(pbmc)<-pbmc$group

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

monocle_cds <-orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 0.75) 
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 0.75)
plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 0.75)
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