# Demo code for manuscript Evaluation of Treatment Effect Modification by Biomarkers Measured Pre- and Post-randomization in the Presence of Non-monotone Missingness
# It first generates a simulated data with similar design as the real dengue data used in the manuscript Fig 3
# Then it illustrates how the code could be used to analyze this simulated data set. 
# Date      Programmer   	
# August12-2020  Y Zhuang      
##########################################################################
library(MASS)
library(mvtnorm)
library(nnet)
library(splines)
rm(list=ls(all=TRUE))

#set your working directory
source("myFunctions.R")
source("SetParameters.R")

#for generating data
myseed<-10
data<-generateData(myseed=myseed,beta_true=beta_true,nuisance_true=nuisance_true)

#get uniform Su across any datasets:
set.seed(1985)
Su<-rnorm(504,mean=3.2,sd=0.75)
Su<-unique(sort(Su[Su>=1&Su<=5])) #length=500
Su[1]<-c

beta.start <- c(1.5, 0.16,  -0.34 , -0.21, -0.25, 0.2360,0.1060,0.2000)
beta.results<-getBeta(data, beta.start=beta.start)

VE.S.marginal<-getVE(data,Su,results=beta.results,type="marginal")
VE.S.Bpos<-getVE(data,Su,results=beta.results,type="Bpos")
VE.S.Bneg<-getVE(data,Su,results=beta.results,type="Bneg")

#Perturbation
VE.per.Mar<-NULL
VE.per.Bpos<-NULL
VE.per.Bneg<-NULL

iter<-500

for (i in 1:iter){
  
  xi<-rexp(nrow(data),rate=1)
  
  beta.results<-getBeta.Perturb(data,beta.start,xi)
  
  VE.S.marginal.per<-getVE.Perturb(data,Su,results.Perturb=beta.results,type="marginal",xi=xi)
  VE.S.Bpos.per<-getVE.Perturb(data,Su,results.Perturb=beta.results,type="Bpos",xi=xi)
  VE.S.Bneg.per<-getVE.Perturb(data,Su,results.Perturb=beta.results,type="Bneg",xi=xi)
  
  VE.per.Mar<-rbind(VE.per.Mar,VE.S.marginal.per)
  VE.per.Bpos<-rbind(VE.per.Bpos,VE.S.Bpos.per)
  VE.per.Bneg<-rbind(VE.per.Bneg,VE.S.Bneg.per)

  rm(VE.S.marginal.per,VE.S.Bpos.per,VE.S.Bneg.per)
}

#Marginal 
VE.S<-VE.S.marginal
VE.per <- VE.per.Mar

SE.RR=apply(log(1-VE.per),2,sd,na.rm=T)
U=rep(NA,nrow(VE.per))
for (kk in 1:length(U)){
  U[kk]=max(abs(log(1-VE.per[kk,])-log(1-VE.S))/SE.RR,na.rm=T)
}
c.alpha=quantile(U,0.95,na.rm=T)

high.ci.Mar=1-exp(log(1-VE.S)-qnorm(0.975)*SE.RR)
low.ci.Mar=1-exp(log(1-VE.S)+qnorm(0.975)*SE.RR)
high.cb.Mar=1-exp(log(1-VE.S)-c.alpha*SE.RR)
low.cb.Mar=1-exp(log(1-VE.S)+c.alpha*SE.RR)

#Bpos
VE.S<-VE.S.Bpos
VE.per <- VE.per.Bpos

SE.RR=apply(log(1-VE.per),2,sd,na.rm=T)
U=rep(NA,nrow(VE.per))
for (kk in 1:length(U)){
  U[kk]=max(abs(log(1-VE.per[kk,])-log(1-VE.S))/SE.RR,na.rm=T)
}
c.alpha=quantile(U,0.95,na.rm=T)

high.ci.Bpos=1-exp(log(1-VE.S)-qnorm(0.975)*SE.RR)
low.ci.Bpos=1-exp(log(1-VE.S)+qnorm(0.975)*SE.RR)
high.cb.Bpos=1-exp(log(1-VE.S)-c.alpha*SE.RR)
low.cb.Bpos=1-exp(log(1-VE.S)+c.alpha*SE.RR)

#Bneg
VE.S<-VE.S.Bneg
VE.per <- VE.per.Bneg

SE.RR=apply(log(1-VE.per),2,sd,na.rm=T)
U=rep(NA,nrow(VE.per))
for (kk in 1:length(U)){
  U[kk]=max(abs(log(1-VE.per[kk,])-log(1-VE.S))/SE.RR,na.rm=T)
}
c.alpha=quantile(U,0.95,na.rm=T)

high.ci.Bneg=1-exp(log(1-VE.S)-qnorm(0.975)*SE.RR)
low.ci.Bneg=1-exp(log(1-VE.S)+qnorm(0.975)*SE.RR)
high.cb.Bneg=1-exp(log(1-VE.S)-c.alpha*SE.RR)
low.cb.Bneg=1-exp(log(1-VE.S)+c.alpha*SE.RR)


save(Su,VE.S.marginal,VE.S.Bpos,VE.S.Bneg,
     high.ci.Mar,low.ci.Mar,high.cb.Mar,low.cb.Mar,
     high.ci.Bpos,low.ci.Bpos,high.cb.Bpos,low.cb.Bpos,
     high.ci.Bneg,low.ci.Bneg,high.cb.Bneg,low.cb.Bneg,
     file=paste0("Results.RData"))
