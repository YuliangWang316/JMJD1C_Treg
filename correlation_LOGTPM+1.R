library(ggcorrplot)
library(ggthemes)
library(Seurat)
library(dplyr)
library(patchwork)
a<-read.csv("c:/Users/xjmik/Downloads/GSE116347_Gene_count_table_human.csv",sep = ",",header = TRUE)
b<-a[!duplicated(a$gene_symbol),]
rownames(b)<-b$gene_symbol
b<-b[,-1]
c<-b[,1:11]
pbmc <- CreateSeuratObject(counts = c, project = "pbmc3k")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc)
d<-pbmc@assays$RNA@data
d<-as.data.frame(d)
corr<-round(cor(t(d)),3)
corrtest<-round(cor_pmat(t(c)),3)
d<-as.data.frame(corr[,"JMJD1C"])
write.table(d,"correlationGSE116347_logTPM+1.txt",sep = "\t")
