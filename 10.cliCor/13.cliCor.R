

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")

install.packages("ggpubr")

library(limma)
library(ggpubr)

setwd("D:\\auto-tou\\13.cliCor")                  
file="Riskscore.txt"                                                     
rt=read.table(file,sep="\t",header=T,check.names=F,row.names=1)            
cli=read.table("N.txt",sep="\t",header=T,check.names=F,row.names=1)    
gene=colnames(rt)[1]
clinical=colnames(cli)[1]


outTab=data.frame()

	rt1=rt
	
	data=cbind(rt1,gene=rt1[,gene])
	data=as.matrix(data[,c(gene,"gene")])
	if(nchar(row.names(data)[1])!=nchar(row.names(cli)[1])){
		row.names(data)=gsub(".$","",row.names(data))}
	data=avereps(data)
	sameSample=intersect(row.names(data),row.names(cli))
	sameData=data[sameSample,]
	sameClinical=cli[sameSample,]
	cliExpData=cbind(as.data.frame(sameClinical),sameData)
	if(nrow(cliExpData)==0){next}
	
	group=levels(factor(cliExpData$sameClinical))
	comp=combn(group,2)
	my_comparisons=list()
    for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
	
	boxplot=ggboxplot(cliExpData, x="sameClinical", y="gene", color="sameClinical",
	          xlab=clinical,
	          ylab=paste(gene,"expression"),
	          legend.title=clinical,
	         
	          add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	pdf(file=paste0(clinical,".",".pdf"),width=5.5,height=5)
	print(boxplot)
	dev.off()



