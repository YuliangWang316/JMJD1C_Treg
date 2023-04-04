library(GEOquery)
gds <- getGEO("GSE14520")
exp<-cbind(gds$`GSE14520-GPL3921_series_matrix.txt.gz`@assayData$exprs,as.data.frame(gds$`GSE14520-GPL3921_series_matrix.txt.gz`@featureData@data$`Gene Symbol`))
colnames(exp)[length(colnames(exp))]<-"Symbol"
gene_symbol<-unique(exp$Symbol)
library(dplyr)
a<-filter(exp,Symbol == gene_symbol[1])
a<-a[,-length(colnames(a))]
b<-as.data.frame(t(as.data.frame(apply(a, MARGIN = 2, median))))
rownames(b)<-gene_symbol[1]
remove(a)
for (i in 2:length(gene_symbol)) {
  a<-filter(exp,Symbol == gene_symbol[i])
  a<-a[,-length(colnames(a))]
  c<-as.data.frame(t(as.data.frame(apply(a, MARGIN = 2, median))))
  remove(a)
  rownames(c)<-gene_symbol[i]
  b<-rbind(b,c)
  remove(c)
}
library(GSVA)
genset<-read.table("c:/Users/xjmik/Downloads/Treg_1C_signature.txt",sep = "\t",header = TRUE)
geneset_list<-as.list(genset)
remove(exp,genset,gene_symbol)

exp<-cbind(gds$`GSE14520-GPL571_series_matrix.txt.gz`@assayData$exprs,as.data.frame(gds$`GSE14520-GPL571_series_matrix.txt.gz`@featureData@data$`Gene Symbol`))
colnames(exp)[length(colnames(exp))]<-"Symbol"
gene_symbol<-unique(exp$Symbol)
library(dplyr)
a<-filter(exp,Symbol == gene_symbol[1])
a<-a[,-length(colnames(a))]
d<-as.data.frame(t(as.data.frame(apply(a, MARGIN = 2, median))))
rownames(d)<-gene_symbol[1]
remove(a)
for (i in 2:length(gene_symbol)) {
  a<-filter(exp,Symbol == gene_symbol[i])
  a<-a[,-length(colnames(a))]
  c<-as.data.frame(t(as.data.frame(apply(a, MARGIN = 2, median))))
  remove(a)
  rownames(c)<-gene_symbol[i]
  d<-rbind(d,c)
  remove(c)
}
remove(gene_symbol,exp)
remove(i)
e<-intersect(rownames(b),rownames(d))
b_new<-b[e,]
d_new<-d[e,]
f<-cbind(b_new,d_new)
remove(b,b_new,d,d_new,e)

gsva<-gsva(expr = as.matrix(f),gset.idx.list = geneset_list,kcdf="Poisson",parallel.sz=20)
gsva<-as.data.frame(t(gsva))
remove(f,geneset_list)
remove(gds)
a<-read.table("c:/Users/xjmik/Downloads/GSE14520_Extra_Supplement.txt",sep = "\t",header = TRUE,row.names = 1)
rownames(a)<-a$Affy_GSM
a<-a[rownames(gsva),]
a<-cbind(a,gsva)
remove(gsva)
b<-a[,c(18,19,22)]
remove(a)
colnames(b)<-c("events","OS","Treg_1C")
library(survival)
library(survminer)
library(dplyr)
df<-b
remove(b)
df<-na.omit(df)
res.cut<-surv_cutpoint(df,time = "OS",event = "events",variables = "Treg_1C" )
summary(res.cut)
res.cat<-surv_categorize(res.cut)
fit<-survfit(Surv(OS,events)~ Treg_1C,data = res.cat)
ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
res.cox<-coxph(Surv(OS,events) ~ Treg_1C,data=res.cat)
summary(res.cox)
test.ph<-cox.zph(res.cox)
ggcoxzph(test.ph)