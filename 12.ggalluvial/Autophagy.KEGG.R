
install.packages("colorspace")
install.packages("stringi")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DOSE")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("clusterProfiler")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("enrichplot")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("C:\\Users\\lexb4\\Desktop\\Autophagy\\13.KEGG")            
rt=read.table("id.txt",sep="\t",header=T,check.names=F)       
rt=rt[is.na(rt[,"entrezID"])==F,]                           
gene=rt$entrezID


kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)  
write.table(kk,file="KEGGId.txt",sep="\t",quote=F,row.names = F)                         


pdf(file="barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = 30)
dev.off()


pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 30)
dev.off()

