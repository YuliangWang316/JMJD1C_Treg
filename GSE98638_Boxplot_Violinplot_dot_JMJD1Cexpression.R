library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(cowplot)
library(tidyverse)
library(beeswarm)
library(ggbeeswarm)
data<-read.table("c:/Users/xjmik/Downloads/Dataset-2023-04-16.txt",sep = "\t",header = TRUE)
data$Group<-factor(data$Group,levels = c("PbmcTr","TumorTr","PbmcTc","TumorTc"))
ggplot(data,aes(x=Group,y=JMJD1C))+ geom_quasirandom(size=0.5,alpha=0.8,aes(color=Group)) + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) +geom_boxplot(alpha=0) + scale_color_manual(values = c("Gray","Gray","Gray","Gray"))
 
# beeswarm(JMJD1C ~ Group, data = data,
#          pch = 16, 
#          xlab = "", ylab = "",
#          cex=0.1)
# bxplot(JMJD1C ~ Group, data = data, add = TRUE)
