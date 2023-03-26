library(dplyr)
library(Seurat)

pbmc.data<-read.table("c:/Users/xjmik/Downloads/GSE99254/GSE99254_NSCLC.TCell.S12346.count.txt",sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
pbmc.data<-pbmc.data[!duplicated(pbmc.data$symbol),]
pbmc.data<-na.omit(pbmc.data)
rownames(pbmc.data)<-pbmc.data$symbol
pbmc.data<-pbmc.data[,-1]
pbmc.metadata<-read.table("c:/Users/xjmik/Downloads/GSE99254/metadata.txt",sep = "\t",header = TRUE,row.names = 1)
pbmc.metadata<-pbmc.metadata[colnames(pbmc.data),]
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata)
remove(pbmc.data,pbmc.metadata)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
remove(all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
Idents(pbmc)<-pbmc@meta.data$majorCluster


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
  if(pbmc@meta.data$sampleType[i] == "TTY"){
    pbmc@meta.data$Bulksample[i] <- "T"
  }
  if(pbmc@meta.data$sampleType[i] == "NTY"){
    pbmc@meta.data$Bulksample[i] <- "N"
  }
  if(pbmc@meta.data$sampleType[i] == "PTY"){
    pbmc@meta.data$Bulksample[i] <- "P"
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
  if(pbmc@meta.data$sampleType[i] == "PTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
  }
  if(pbmc@meta.data$sampleType[i] == "TTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
  }
  if(pbmc@meta.data$sampleType[i] == "NTY"){
    pbmc@meta.data$Bulksample2[i] <- "TY"
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


Idents(pbmc)<-pbmc@meta.data$majorCluster

Tregtumor<-subset(pbmc,idents = c("CD4_C9-CTLA4"))
Tregtumor <- NormalizeData(Tregtumor)
Tregtumor <- FindVariableFeatures(Tregtumor, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tregtumor)
Tregtumor <- ScaleData(Tregtumor, features = all.genes)
remove(all.genes)
Tregtumor <- RunPCA(Tregtumor, features = VariableFeatures(object = Tregtumor),)
Tregtumor <- RunUMAP(Tregtumor,dims = 1:30)
Tregtumor <- RunTSNE(Tregtumor,dims = 1:30)
i<-data.frame(0,0,0,0,0,0,0,0,0,0,0,0)
colnames(i)<-c("a","b","c","d","e","f","g","h","j","k","l","m")
for (a in seq(2,50,1)) {
  for (b in seq(2,100,1)) {
    for (c in seq(10,500,10)) {
      for (d in c("pca")) {
        for (e in c("euclidean","cosine","manhattan","hamming")) {
          for (f in seq(0.1,3,0.1)) {
            for (g in c(1,2,3)) {
              for (h in 1) {
                temp    <- try(Tregtumor_new <- FindNeighbors(Tregtumor, dims = 1:a,k.param = 20,n.trees = c,reduction = d,annoy.metric = e),silent=FALSE)
                if(!('try-error' %in% class(temp))){
                  Tregtumor_new <- temp
                  temp    <- try(Tregtumor_new <- FindClusters(Tregtumor_new, resolution = f,algorithm = g,modularity.fxn = h ),silent=FALSE)
                  if(!('try-error' %in% class(temp))){
                    Tregtumor_new <- temp
                    if(length(unique(Tregtumor_new$seurat_clusters)) < 10){
                      temp    <- try(Tregtumor.marker<-FindAllMarkers(Tregtumor_new,only.pos = TRUE),silent=FALSE)
                      if(!('try-error' %in% class(temp))){  
                        Tregtumor.marker<-temp
                        for (n in 1:length(rownames(Tregtumor.marker))) {
                          if(Tregtumor.marker[n,7] %in% "JMJD1C"){
                            if(Tregtumor.marker[n,5] < 0.05){
                              k<-data.frame(a,b,c,d,e,f,g,h,Tregtumor.marker[n,1],Tregtumor.marker[n,2],Tregtumor.marker[n,5],Tregtumor.marker[n,6])
                              colnames(k)<-c("a","b","c","d","e","f","g","h","j","k","l","m")
                              i<-rbind(i,k)
                            }
                          }
                        }
                      }else{
                        next
                      }
                    }
                  }else{
                    next
                  }
                }else{
                  next
                }
              }
            }
          }
        }
      }
    }
  }
}