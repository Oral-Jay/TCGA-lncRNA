
setwd("C:\\Users\\4_Calibration")

library(rms)
library(foreign)
library(survival)


tcga<-read.table("clinical.txt",header=T,sep="\t")
 
tcga$age<-factor(tcga$age,labels=c("<50","50-59","60-69",">=70"))
tcga$sex<-factor(tcga$sex,labels=c("FEMALE","MALE"))
tcga$race<-factor(tcga$race,labels=c("WHITE","BLACK OR AFRICAN AMERICAN","ASIAN","AMERICAN INDIAN OR ALASKA NATIVE"))
tcga$smoking<-factor(tcga$smoking,labels=c("1","2","3","4","5"))
tcga$radiation<-factor(tcga$radiation,labels=c("YES","NO"))
tcga$pharmaceutical<-factor(tcga$pharmaceutical,labels=c("YES","NO"))
tcga$stage_T<-factor(tcga$stage_T,labels=c("T1","T1a","T1b","T1c","T2","T2a","T2b","T3","T4","TX"))
tcga$stage_M<-factor(tcga$stage_M,labels=c("M0","M1","M1a","M1b","MX"))
tcga$stage_N<-factor(tcga$stage_N,labels=c("N0","N1","N2","N3","NX"))
tcga$surgery<-factor(tcga$surgery,labels=c("YES","NO"))

cox1 <- cph(Surv(futime,fustat)~age + gender + grade + stage + T + N + riskScore,surv=T,x=T, y=T,time.inc = 1*365*3,data=tcga) 

cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*365*3, m= 30, B=1000)


pdf("calibrate3.pdf",12,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1),ylim = c(0,1),xlab ="Nomogram Predicted Survival",ylab="Actual Survival",col="blue")
lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 3,col ="black" ,pch = 16)
mtext('')
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()

