pbmc<-readRDS("d:/PanT.rds")
Treg_T<-readRDS("d:/PanTreg_tumor.rds")
library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
JMJD1Chi<-subset(Treg_T,ident = c("1","5","7"))
name<-colnames(JMJD1Chi)
remove(JMJD1Chi,Treg_T)
for (i in 1:length(name)) {
  for (j in 1:length(rownames(pbmc@meta.data))) {
    if(name[i] == colnames(pbmc)[j]){
      pbmc@meta.data$Group[j] == "hi"
    }else{
      pbmc@meta.data$Group[j] == "lo"
    }
  }
}
