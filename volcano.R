setwd("D:/")
diff_gene <- read.table("volcanoplot.txt", sep="\t", header=TRUE, row.names=1)
diff_gene=as.data.frame(diff_gene)
gene_list=diff_gene[,c("log2FoldChange","padj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 1 & gene_list$padj < 0.05)
colored_point<-gene_list[gene_list$threshold == "TRUE",]

gene_list$threshold<-as.character(gene_list$threshold)

Mycolors<-c("Black","Black")
library("ggplot2")
pdf("vocano.pdf")

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=threshold)) + geom_point(alpha=0.4, size=1.75)  + xlim(c(-3, 3)) + ylim(c(0, 3)) +xlab("log2 fold change") + ylab("-log10 p-value")+ theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) + scale_color_manual(values = Mycolors)
print(g)
dev.off()
