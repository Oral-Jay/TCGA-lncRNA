
install.packages('survival')

library(survival)                                                  
setwd("D:\\auto-tou\\12.multicox")           
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)    

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)


outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)


riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
    file="risk.txt",
    sep="\t",
    quote=F,
    row.names=F)

