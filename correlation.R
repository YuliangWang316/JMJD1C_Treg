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
corr<-round(cor(t(c)),3)
corrtest<-round(cor_pmat(t(c)),3)
d<-as.data.frame(corr[,"JMJD1C"])
write.table(d,"GSE116347_Treg_1c_Correlation.txt",sep = "\t")
