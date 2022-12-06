
setwd("D:\\auto-tou\\1-3years")

library(rms)
library(foreign)
library(survival)


tcga<-read.table("clinical.txt",header=T,sep="\t")
names(tcga)


#tcga$age<-factor(tcga$age,labels=c("<50","50-59","60-69",">=70"))
tcga$gender<-factor(tcga$gender,levels=c(0,1),labels=c("FEMALE","MALE"))
tcga$grade<-factor(tcga$grade,labels=c("1","2","3",'4'))
#tcga$smoking<-factor(tcga$smoking,labels=c("1","2","3","4","5"))
#tcga$radiation<-factor(tcga$radiation,labels=c("YES","NO"))
#tcga$pharmaceutical<-factor(tcga$pharmaceutical,labels=c("YES","NO"))
tcga$stage<-factor(tcga$stage,levels=c(1,2,3,4),labels=c("1","2","3",'4'))
tcga$T<-factor(tcga$T,labels=c("1","2","3",'4'))
tcga$N<-factor(tcga$N,labels=c("N0","N1","N2"))
#tcga$surgery<-factor(tcga$surgery,labels=c("YES","NO"))


ddist <- datadist(tcga)
options(datadist='ddist')


cox <- cph(Surv(futime,fustat) ~age + gender + grade + stage + T + N + riskScore,surv=T,x=T, y=T,data=tcga) 
surv <- Survival(cox)

surv <- Survival(cox)
sur_3_year<-function(x)surv(1*365*3,lp=x)
sur_1_year<-function(x)surv(1*365*1,lp=x)
nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year),lp= F,funlabel=c('1-Year Survival','3-Year survival'),maxscale=100,fun.at=c('0.9','0.7','0.5','0.3',"0.1"))


pdf("nom.pdf",15,10)
plot(nom_sur,xfrac=0.25)
dev.off()
nom_sur
