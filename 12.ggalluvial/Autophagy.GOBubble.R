
install.packages("digest")
install.packages("GOplot")

library(GOplot)
setwd("G:\\BaiduNetdiskDownload\\12.GOBubble")              

ego=read.table("GO.txt", header = T,sep="\t",check.names=F)          
go=data.frame(Category = ego$ONTOLOGY,ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)


id.fc <- read.table("id.txt", header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$gene, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]
circ <- circle_dat(go, genelist)


pdf(file="GOBubble.pdf",width = 10,height = 8)
GOBubble(circ, labels = 3,table.legend =F)
dev.off()


pdf(file="GOCircle.pdf",width = 11,height = 6)
GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10)          
dev.off()


termNum = 20                                     
geneNum = nrow(genelist)                         
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="GOHeat.pdf",width = 11,height = 5)
GOHeat(chord, nlfc =1, fill.col = c('red', 'white', 'blue'))
dev.off()
