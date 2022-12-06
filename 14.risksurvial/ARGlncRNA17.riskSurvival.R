
install.packages("survival")

library(survival)
setwd("D:\\auto-tou\\17.risksurvial")       

data=read.table("risk.txt",header=T,sep="\t",check.names=F)
diff=survdiff(Surv(futime, fustat) ~risk,data = data)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = data)
pdf(file="survival.pdf",width=5.5,height=5)
plot(fit,
	 lwd=2,
	 col=c("red","blue"),
	 xlab="Time (year)",
	 ylab="Survival rate",
	 main=paste("Survival curve (p=", pValue ,")",sep=""),
	 mark.time=T)
legend("topright",
	 c("High risk", "Low risk"),
     lwd=2,
	 col=c("red","blue"))
dev.off()

summary(fit)


