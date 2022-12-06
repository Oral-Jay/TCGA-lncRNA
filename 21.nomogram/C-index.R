
install.packages("rms") 
install.packages("foreign") 
install.packages("survival") 


setwd("C:\\Users\\scikuangren\\Desktop\\Clinical_predictive_model\\3_C-index")

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


ddist <- datadist(tcga)
options(datadist='ddist')


fmla1 <- as.formula(Surv(futime,fustat) ~age + gender + grade + stage + T + N + riskScore)
cox2 <- cph(fmla1,data=tcga)

coxpe <- predict(cox2)
c_index=1-rcorr.cens(coxpe,Surv(tcga$futime,tcga$fustat))
c_index

