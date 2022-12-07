#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)


setwd("D:\\auto-tou\\08.ARGexp")


rt=read.table("protein.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]


gene=read.table("gene.txt", header=F, check.names=F, sep="\t")
sameGene=intersect(as.vector(gene[,1]),rownames(data))
geneExp=data[sameGene,]


out=rbind(ID=colnames(geneExp),geneExp)
write.table(out,file="ARGexp.txt",sep="\t",quote=F,col.names=F)

