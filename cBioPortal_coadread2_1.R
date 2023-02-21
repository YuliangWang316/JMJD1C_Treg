library(cBioPortalData)
library(AnVIL)
library(GSVA)
library(TCGAbiolinks)
lihc_genexpression<-read.table("c:/Users/xjmik/Downloads/coadread_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt",sep = "\t",header = TRUE,check.names = FALSE)
lihc_genexpression<-lihc_genexpression[!duplicated(lihc_genexpression$Hugo_Symbol),]
rownames(lihc_genexpression)<-lihc_genexpression[,1]
lihc_genexpression<-lihc_genexpression[,-1]
lihc_genexpression<-lihc_genexpression[,-1]
lihc_genexpression<-na.omit(lihc_genexpression)
lihc_clinical<-read.table("c:/Users/xjmik/Downloads/coadread_tcga_pan_can_atlas_2018/clinical_combine.txt",sep = "\t",header = TRUE,row.names = 1)
samplesTP <- TCGAquery_SampleTypes(colnames(lihc_genexpression), typesample = c("TP"))
lihc_genexpression<-lihc_genexpression[,samplesTP]
genset<-read.table("c:/Users/xjmik/Downloads/1C_signature.txt",sep = "\t")
colnames(genset)<-"J1C_signature"
geneset_list<-as.list(genset)

    gsva_gl<-gsva(expr = as.matrix(lihc_genexpression),gset.idx.list = geneset_list,kcdf="Gaussian",parallel.sz=20)
    gsva_gl<-as.data.frame(t(gsva_gl))
    gsva_gl$sample <- sapply(strsplit(rownames(gsva_gl),'-'),function(x) paste0(x[1:3],collapse="-"))
    
    
    lihc_clinical<-lihc_clinical[gsva_gl$sample,]
    lihc_clinical$"J1C_signature" <- gsva_gl[match(rownames(lihc_clinical),gsva_gl$sample),][,"J1C_signature"]
    
    
    df<-subset(lihc_clinical,select =c("OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS","J1C_signature","PFS_STATUS","PFS_MONTHS","DSS_STATUS","DSS_MONTHS"))
    
    df <- df[!is.na(df$J1C_signature),]
    
    df<-df[which(df$OS_MONTHS!="[Not Available]"),]
    df<-df[which(df$OS_MONTHS!= 0),]
    
    
    for (j in 1:length(rownames(df))) {
      if(df$OS_STATUS[j] == "0:LIVING"){
        df$events[j]<-0
      }else if(df$OS_STATUS[j] == "1:DECEASED"){
        df$events[j]<-1
      }
    }
    library(survival)
    library(survminer)
    df$OS_MONTHS<-as.numeric(df$OS_MONTHS)
    res.cut<-surv_cutpoint(df,time = "OS_MONTHS",event = "events",variables = c("J1C_signature"))
    summary(res.cut)
    res.cat<-surv_categorize(res.cut)
    fit<-survfit(Surv(OS_MONTHS,events)~ J1C_signature ,data = res.cat)
    ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
    
    df<-subset(lihc_clinical,select =c("OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS","J1C_signature","PFS_STATUS","PFS_MONTHS","DSS_STATUS","DSS_MONTHS"))
    
    df <- df[!is.na(df$J1C_signature),]
    
    df<-df[which(df$DFS_MONTHS!="[Not Available]"),]
    df<-df[which(df$DFS_MONTHS!= 0),]
    
    
    for (j in 1:length(rownames(df))) {
      if(df$DFS_STATUS[j] == "0:DiseaseFree"){
        df$events[j]<-0
      }else if(df$DFS_STATUS[j] == "1:Recurred/Progressed"){
        df$events[j]<-1
      }
    }
    library(survival)
    library(survminer)
    df$DFS_MONTHS<-as.numeric(df$DFS_MONTHS)
    res.cut<-surv_cutpoint(df,time = "DFS_MONTHS",event = "events",variables = c("J1C_signature"))
    summary(res.cut)
    res.cat<-surv_categorize(res.cut)
    fit<-survfit(Surv(DFS_MONTHS,events)~ J1C_signature ,data = res.cat)
    ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
    
    df<-subset(lihc_clinical,select =c("OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS","J1C_signature","PFS_STATUS","PFS_MONTHS","DSS_STATUS","DSS_MONTHS"))
    
    df <- df[!is.na(df$J1C_signature),]
    
    df<-df[which(df$PFS_MONTHS!="[Not Available]"),]
    df<-df[which(df$PFS_MONTHS!= 0),]
    
    
    for (j in 1:length(rownames(df))) {
      if(df$PFS_STATUS[j] == "0:CENSORED"){
        df$events[j]<-0
      }else if(df$PFS_STATUS[j] == "1:PROGRESSION"){
        df$events[j]<-1
      }
    }
    library(survival)
    library(survminer)
    df$PFS_MONTHS<-as.numeric(df$PFS_MONTHS)
    res.cut<-surv_cutpoint(df,time = "PFS_MONTHS",event = "events",variables = c("J1C_signature"))
    summary(res.cut)
    res.cat<-surv_categorize(res.cut)
    fit<-survfit(Surv(PFS_MONTHS,events)~ J1C_signature  ,data = res.cat)
    ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
    
    df<-subset(lihc_clinical,select =c("OS_STATUS","OS_MONTHS","DFS_STATUS","DFS_MONTHS","J1C_signature","PFS_STATUS","PFS_MONTHS","DSS_STATUS","DSS_MONTHS"))
    
    df <- df[!is.na(df$J1C_signature),]
    
    df<-df[which(df$DSS_MONTHS!="[Not Available]"),]
    df<-df[which(df$DSS_MONTHS!= 0),]
    df<-df[!is.na(df$DSS_STATUS),]
    
    for (j in 1:length(rownames(df))) {
      if(df$DSS_STATUS[j] == "0:ALIVE OR DEAD TUMOR FREE"){
        df$events[j]<-0
      }else if(df$DSS_STATUS[j] == "1:DEAD WITH TUMOR"){
        df$events[j]<-1
      }
    }
    library(survival)
    library(survminer)
    df$DSS_MONTHS<-as.numeric(df$DSS_MONTHS)
    res.cut<-surv_cutpoint(df,time = "DSS_MONTHS",event = "events",variables = c("J1C_signature"))
    summary(res.cut)
    res.cat<-surv_categorize(res.cut)
    fit<-survfit(Surv(DSS_MONTHS,events)~ J1C_signature  ,data = res.cat)
    ggsurvplot(fit,data = res.cat,risk.table = TRUE,pval = TRUE)
    

