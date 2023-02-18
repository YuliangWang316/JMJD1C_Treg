library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(EDASeq)
library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(org.Hs.eg.db)
library(tidyverse)
library(GSVA)
library(IDConverter)
library("limma")
setwd("D:/TEST3/31.CIBERSORT/")
request_cancer=c("ACC","CHOL","BLCA","COAD","READ","BRCA","LGG","GBM","PCPG","CESC","ESCA","STAD","UVM","HNSC","KIRC","KIRP","KICH","LIHC","LUAD","LUSC","DLBC","LAML","OV","PAAD","MESO","PRAD","SKCM","SARC","TGCT","THYM","THCA","UCEC","UCS")
i<-"LIHC"
cancer_type=paste("TCGA",i,sep="-")
print(cancer_type)
clinical<-GDCquery_clinic(project = cancer_type, type = "clinical")
query <- GDCquery(project = cancer_type, 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

GDCdownload(query, method = "api", files.per.chunk = 100)
expdat <- GDCprepare(query = query)
count_matrix<-as.data.frame(assay(expdat))
count_gl<-TCGAanalyze_Normalization(count_matrix, geneInfoHT,method =  'geneLength')
remove(count_matrix,expdat,query)
genename<-rownames(count_gl)
e<-bitr(geneID = genename,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
e<-e[!duplicated(e$SYMBOL),]
count_gl<-count_gl[e$ENSEMBL,]
rownames(count_gl)<-e$SYMBOL
remove(e,genename)
count_gl<-as.data.frame(count_gl)

newname_gl<-filter_tcga_barcodes(colnames(count_gl),analyte_target = "RNA")
count_gl_new<-count_gl[,newname_gl]
count_gl_new_new<-count_gl_new[rowMeans(count_gl_new)>0,]


library(tidyverse)
library(DESeq2)

samplesTP <- TCGAquery_SampleTypes(colnames(count_gl_new_new), typesample = c("TP"))
count_gl_new_new_tumor<-count_gl_new_new[,samplesTP]
for (p in 1:length(colnames(count_gl_new_new_tumor))) {
  for (q in 1:length(rownames(count_gl_new_new_tumor))) {
    sum_gene<-sum(count_gl_new_new_tumor[,p])
    count_gl_new_new_tumor[q,p]<-count_gl_new_new_tumor[q,p] / sum_gene * 10E10
  }
}
genset<-read.table("c:/Users/xjmik/Downloads/1C_signature.txt",sep = "\t")
colnames(genset)<-"J1C_signature"
geneset_list<-as.list(genset)
gsva_gl<-gsva(expr = as.matrix(count_gl_new_new_tumor),gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20)
gsva_gl<-as.data.frame(t(gsva_gl))
gsva_gl$sample <- sapply(strsplit(rownames(gsva_gl),'-'),function(x) paste0(x[1:3],collapse="-"))
newname_gl<-filter_tcga_barcodes(rownames(gsva_gl),analyte_target = "RNA")
remove(count_gl,count_gl_new,geneset_list,genset,samplesTP)
gsva_gl<-gsva_gl[newname_gl,]
remove(newname_gl)
remove(count_gl_new_new)
clinical$"Treg_1C_gl" <- gsva_gl[match(clinical$submitter_id,gsva_gl$sample),][,"J1C_signature"]

v <-voom(count_gl_new_new_tumor, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol_LIHC.txt",sep="\t",quote=F,col.names=F)  
source("TMEimmune31.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol_LIHC.txt", perm=100, QN=TRUE)
cibersort_tumor<-results
remove(results)
cibersort_tumor1<-data.frame(rownames(cibersort_tumor),cibersort_tumor[,9])
rownames(cibersort_tumor1)<-rownames(cibersort_tumor)
colnames(cibersort_tumor1)<-c("name","Treg")
cibersort_tumor<-cibersort_tumor1
remove(cibersort_tumor1)
remove(out,v,i,p,q,sum_gene,request_cancer)
cibersort_tumor$sample <- sapply(strsplit(rownames(cibersort_tumor),'-'),function(x) paste0(x[1:3],collapse="-"))
K <- cibersort_tumor[match(clinical$submitter_id,cibersort_tumor$sample),]
clinical<-cbind(clinical,K)
remove(K)
# library(estimate)
# write.table(count_gl_new_new_tumor,file="uniq.symbol.txt",sep="\t",quote=F,col.names=TRUE)
# filterCommonGenes(input.f="uniq.symbol.txt", 
#                   output.f="commonGenes.gct", 
#                   id="GeneSymbol")
# 
# estimateScore(input.ds = "commonGenes.gct",
#               output.ds="estimateScore.gct", 
#               platform="illumina")
# scores=read.table("estimateScore.gct",skip = 2,header = T)
# rownames(scores)=scores[,1]
# scores=t(scores[,3:ncol(scores)])
# rownames(scores)=gsub("\\.","\\-",rownames(scores))
# out=rbind(ID=colnames(scores),scores)
# remove(scores)
# out<-out[-1,]
remove(cibersort_tumor,gsva_gl)
df<-subset(clinical,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,Treg_1C_gl,age_at_index,gender,ajcc_pathologic_stage,Treg))
df <- df[!is.na(df$Treg_1C_gl),]
df <- df[!is.na(df$Treg),]

df<-df[which(df$vital_status!="NA"),]
for (j in 1:length(rownames(df))) {
  if(is.na(df$days_to_death[j])){
    df$Time[j] <- df$days_to_last_follow_up[j]
  }else if(is.na(df$days_to_last_follow_up[j]) ){
    df$Time[j] <- df$days_to_death[j]
  }
  else if(df$days_to_death[j] >=df$days_to_last_follow_up[j]){
    df$Time[j] <-df$days_to_death[j]
  }
}
df<-df[which(df$Time != 0),]
for (j in 1:length(rownames(df))) {
  if(df$vital_status[j] == "Alive"){
    df$events[j]<-0
  }else if(df$vital_status[j] == "Dead"){
    df$events[j]<-1
  }
}
res.cut<-surv_cutpoint(df,time = "Time",event = "events",variables = c("Treg_1C_gl","Treg" ))
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(Time,events)~ Treg_1C_gl + Treg ,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
res.cat_new<-res.cat[which(res.cat$Treg == "low"),]
fit_new<-survfit(Surv(Time,events)~ Treg_1C_gl ,data = res.cat_new)
ggsurvplot(fit_new,data = res.cat_new,risk.table = TRUE,pval = TRUE)
