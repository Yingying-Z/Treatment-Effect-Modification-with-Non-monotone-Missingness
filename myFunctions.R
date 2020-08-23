#########################################################################
# DESCRIPTION: functions that will be called in Dissertation_2 simulation
#
# CODED BY: Yingying Zhuang
# SAVED AS: myFunctions.R
#
# HISTORY:
# Date      Programmer   	
# Mar18-2016  Y Zhuang      Version 1: 
# Apr26-2016  Y Zhuang      Version 2: 
##########################################################################
generateData <- function(myseed,beta_true,nuisance_true){
  alpha0<-nuisance_true[[1]]
  alpha1<-nuisance_true[[2]]
  sigma1<-nuisance_true[[3]]
  gamma0<-nuisance_true[[4]]
  gamma1<-nuisance_true[[5]]
  gamma2<-nuisance_true[[6]]
  sigma2<-nuisance_true[[7]]
  N<-10000
  W<-rep(c(1,2,3,4),each=N/4)
  data<-data.frame(W)
  data$W2 <- data$W==2
  data$W3 <- data$W==3
  data$W4 <- data$W==4
  
  n1<-round(N/4*(2/3),0)
  data$Z<-rep(c(rep(0,N/4-n1),rep(1,n1)),4)
  
  ######Nuisance parameters
  #nuisance<-list(alpha0,alpha1, sigma1, gamma0, gamma1, gamma2, sigma2, muB, muS, sigmaB, sigmaS, rho)
  set.seed(myseed)
  B1<-rnorm(N/4,mean=alpha0+alpha1[1],sd=sigma1)
  set.seed(myseed)
  B2<-rnorm(N/4,mean=alpha0+alpha1[2],sd=sigma1)
  set.seed(myseed)
  B3<-rnorm(N/4,mean=alpha0+alpha1[3],sd=sigma1)
  set.seed(myseed)
  B4<-rnorm(N/4,mean=alpha0+alpha1[4],sd=sigma1)
  data$B<-c(B1,B2,B3,B4)
  
  mu<-gamma0+gamma1*data$B+gamma2[2]*data$W2+gamma2[3]*data$W3+gamma2[4]*data$W4
  set.seed(myseed)
  data$S<-rnorm(N,mean=mu,sd=sigma2)

  data$S <- ifelse(data$S<1,1,data$S) 
  data$B <- ifelse(data$B<1,1,data$B) 
  
  X<-cbind(rep(1,N),data$Z,data$S,data$S*data$Z,data$B,data$W2,data$W3,data$W4)
  R<-X%*%beta_true
  R<-pnorm(R)
  mean(R[data$Z==0]) #should be about 0.04
  mean(R[data$Z==1]) #should be about 0.02
  
  
  set.seed(myseed)
  data$Y<-rbinom(N,1,R)
  rm(X,R)
  
  #prob. of being sampled for Immunogenicity Set
  prob.FASI<-0.35
  set.seed(myseed+1000)
  data$FASI<-rbinom(N,1,prob.FASI)
  
  #all subjects in the Immunogenicity Set has Baseline titer
  data$B<-ifelse(data$FASI==1,data$B,NA)
  data$delta_B<-data$FASI
  #all subjects in the Immunogenicity Set and all cases has S
  data$S<-ifelse(data$FASI==1|data$Y==1,data$S,NA)
  data$delta_S<-data$FASI==1|data$Y==1
  rm(prob.FASI)
  ################Done with creating simulation data
  return(data)
}



#Estimate B|W

## (Negative) log-likelihood function
#theta=(alpha0,alpha12,alpha13,alpha14,sigma_1)
logl.BonW <- function(theta,data) {
  dat<-data[data$delta_B==1,]
  alpha0<-theta[1]
  alpha1<-c(0,theta[2:4])
  sigma1<-theta[5]
  mu1<-ifelse(dat$W==1,alpha0+alpha1[1],ifelse(dat$W==2,alpha0+alpha1[2],ifelse(dat$W==3,alpha0+alpha1[3],alpha0+alpha1[4])))
  ff<-ifelse(dat$B==c,pnorm(c,mean=mu1,sd=sigma1), dnorm(dat$B,mean=mu1,sd=sigma1))
  return(-sum(log(ff)))
}


#Estimate S|B,W
#theta=c(gamma0,gamma1,gamma22,gamma23,gamma24,sigma2)
logl.SonBW<-function(theta,data,alpha0,alpha1,sigma1){
  dat<-data[data$delta_B==1&data$Z==1,]
  gamma0<-theta[1]
  gamma1<-theta[2]
  gamma2<-c(0,theta[3:5])
  sigma2<-theta[6]
  
  ans<-0
  #S=c,B=c
  dd<-dat[dat$S==c&dat$B==c,] #6
  if(nrow(dd)>0){
  myid<-1:nrow(dd)
  integrand<-function(b,id){
    if(dd$W[id]==1) {mygamma2<-gamma2[1]; mu1<-alpha0+alpha1[1]}
    if(dd$W[id]==2) {mygamma2<-gamma2[2]; mu1<-alpha0+alpha1[2]}
    if(dd$W[id]==3) {mygamma2<-gamma2[3]; mu1<-alpha0+alpha1[3]}
    if(dd$W[id]==4) {mygamma2<-gamma2[4]; mu1<-alpha0+alpha1[4]}
    
    out<-pnorm(c, mean= gamma0+gamma1*b+mygamma2,sd=sigma2)*dnorm(b, mean=mu1,sd=sigma1)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=-Inf,upper=c,id=id)$value
  area2<-Vectorize(area)
  num<-area2(myid)
  
  mu1<-ifelse(dd$W==1,alpha0+alpha1[1],ifelse(dd$W==2,alpha0+alpha1[2],ifelse(dd$W==3,alpha0+alpha1[3],alpha0+alpha1[4])))
  den<-pnorm(c,mean=mu1,sd=sigma1)
  ff<-num/den
  ans<-ans+sum(log(ff))
  }
  
  #S=c,B>c
  dd<-dat[dat$S==c&dat$B>c,] #0
  if(nrow(dd)>0){
  mu2<-ifelse(dd$W==1,gamma0+gamma1*dd$B+gamma2[1],ifelse(dd$W==2,gamma0+gamma1*dd$B+gamma2[2],ifelse(dd$W==3,gamma0+gamma1*dd$B+gamma2[3],gamma0+gamma1*dd$B+gamma2[4])))
  ff<-dnorm(c,mean=mu2,sd=sigma2)
  ans<-ans+sum(log(ff))
  }
  
  #S>c,B=c
  dd<-dat[dat$S>c&dat$B==c,] #101
  if(nrow(dd)>0){
  myid<-1:nrow(dd)
  integrand<-function(b,id){
    if(dd$W[id]==1) {mygamma2<-gamma2[1]; mu1<-alpha0+alpha1[1]}
    if(dd$W[id]==2) {mygamma2<-gamma2[2]; mu1<-alpha0+alpha1[2]}
    if(dd$W[id]==3) {mygamma2<-gamma2[3]; mu1<-alpha0+alpha1[3]}
    if(dd$W[id]==4) {mygamma2<-gamma2[4]; mu1<-alpha0+alpha1[4]}
    
    out<-dnorm(dd$S[id], mean= gamma0+gamma1*b+mygamma2,sd=sigma2)*dnorm(b, mean=mu1,sd=sigma1)
    return(out)
  }
  area<-function(id) integrate(integrand,lower=-Inf,upper=c,id=id)$value
  area2<-Vectorize(area)
  num<-area2(myid)
  
  mu1<-ifelse(dd$W==1,alpha0+alpha1[1],ifelse(dd$W==2,alpha0+alpha1[2],ifelse(dd$W==3,alpha0+alpha1[3],alpha0+alpha1[4])))
  den<-pnorm(c,mean=mu1,sd=sigma1)
  ff<-num/den
  ans<-ans+sum(log(ff))
  }
  
  #S>c,B>c
  dd<-dat[dat$S>c&dat$B>c,] #2304
  mu2<-ifelse(dd$W==1,gamma0+gamma1*dd$B+gamma2[1],ifelse(dd$W==2,gamma0+gamma1*dd$B+gamma2[2],ifelse(dd$W==3,gamma0+gamma1*dd$B+gamma2[3],gamma0+gamma1*dd$B+gamma2[4])))
  ff<-dnorm(dd$S,mean=mu2,sd=sigma2)
  ans<-ans+sum(log(ff))
  
  return(-ans)
}





###Calculate observed likelihood

like<-function(beta,data,nuisance){
  beta_0<-beta[1];beta_1<-beta[2];beta_2<-beta[3];beta_3<-beta[4];beta_4<-beta[5];beta_52<-beta[6];beta_53<-beta[7];beta_54<-beta[8]
  
  alpha0<-nuisance[[1]]
  alpha1<-nuisance[[2]]
  sigma1<-nuisance[[3]]
  gamma0<-nuisance[[4]]
  gamma1<-nuisance[[5]]
  gamma2<-nuisance[[6]]
  sigma2<-nuisance[[7]]
  muB<-nuisance[[8]]
  muS<-nuisance[[9]]
  sigmaB<-nuisance[[10]]
  sigmaS<-nuisance[[11]]
  rho<-nuisance[[12]]
  sigma3<-sigmaB*sqrt(1-rho^2)
  sigma<-matrix(c(sigmaB^2,rho*sigmaB*sigmaS,rho*sigmaB*sigmaS,sigmaS^2),nrow=2)
  
  ans <- 0
  
  # delta_S=delta_B=1,Z=1
  dat<-data[data$delta_S==1&data$delta_B==1&data$Z==1,] #2415
  ff<-ifelse(dat$W==1,beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*dat$B,
             ifelse(dat$W==2,beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*dat$B+beta_52,
                    ifelse(dat$W==3, beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*dat$B+beta_53, beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*dat$B+beta_54)))
  
  risk<-pnorm(ff)
  rr<-ifelse(dat$Y==1,risk,1-risk)
  ans<-ans+sum(log(rr))
  
  #delta_S=delta_B=1,Z=0
  #if B>c 
  dat<-data[data$delta_S==1&data$delta_B==1&data$Z==0&data$B>c,] #1053
  mu2<- ifelse(dat$W==1, gamma0+gamma1*dat$B+gamma2[1],
               ifelse(dat$W==2, gamma0+gamma1*dat$B+gamma2[2],
                      ifelse(dat$W==3, gamma0+gamma1*dat$B+gamma2[3],gamma0+gamma1*dat$B+gamma2[4])))
  ff<-ifelse(dat$W==1,beta_0+beta_2*c+beta_4*dat$B,
             ifelse(dat$W==2,beta_0+beta_2*c+beta_4*dat$B+beta_52,
                    ifelse(dat$W==3, beta_0+beta_2*c+beta_4*dat$B+beta_53, beta_0+beta_2*c+beta_4*dat$B+beta_54)))
  term1<- pnorm( (c-mu2)/sigma2 ) * pnorm(ff)
  
  myid<-1:nrow(dat)
  integrand<-function(s,id){
    if(dat$W[id]==1) {mu2<-gamma0+gamma1*dat$B[id]+gamma2[1]; beta5<-0}
    if(dat$W[id]==2) {mu2<-gamma0+gamma1*dat$B[id]+gamma2[2]; beta5<-beta_52}
    if(dat$W[id]==3) {mu2<-gamma0+gamma1*dat$B[id]+gamma2[3]; beta5<-beta_53}
    if(dat$W[id]==4) {mu2<-gamma0+gamma1*dat$B[id]+gamma2[4]; beta5<-beta_54}
    
    out<-pnorm(beta_0+beta_2*s+beta_4*dat$B[id]+beta5)*dnorm(s, mean=mu2, sd=sigma2)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
  area2<-Vectorize(area)
  
  term2<-area2(myid)
  
  risk<-term1+term2
  rr<-ifelse(dat$Y==1,risk,1-risk)
  ans<-ans+sum(log(rr))
  
  #if B=c
  dat<-data[data$delta_S==1&data$delta_B==1&data$Z==0&data$B==c,] #64
  if(nrow(dat)>0){
  p1<-as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[1],muS[1]),sigma = sigma))/pnorm( (c-muB[1])/sigmaB)
  p1<-c(p1,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[2],muS[2]),sigma = sigma))/pnorm( (c-muB[2])/sigmaB))
  p1<-c(p1,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[3],muS[3]),sigma = sigma))/pnorm( (c-muB[3])/sigmaB))
  p1<-c(p1,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[4],muS[4]),sigma = sigma))/pnorm( (c-muB[4])/sigmaB))
  
  term1<-ifelse(dat$W==1, p1[1]*pnorm(beta_0+beta_2*c+beta_4*c),
                ifelse(dat$W==2, p1[2]*pnorm(beta_0+beta_2*c+beta_4*c+beta_52),
                       ifelse(dat$W==3,p1[3]*pnorm(beta_0+beta_2*c+beta_4*c+beta_53),p1[4]*pnorm(beta_0+beta_2*c+beta_4*c+beta_54))))
  
  myid<-1:nrow(dat)
  integrand<-function(s,id){
    if(dat$W[id]==1) {my.muS<-muS[1]; my.muB<-muB[1]; beta_5<-0}
    if(dat$W[id]==2) {my.muS<-muS[2]; my.muB<-muB[2]; beta_5<-beta_52}
    if(dat$W[id]==3) {my.muS<-muS[3]; my.muB<-muB[3]; beta_5<-beta_53}
    if(dat$W[id]==4) {my.muS<-muS[4]; my.muB<-muB[4]; beta_5<-beta_54}
    
    mu3<-my.muB+rho*(sigmaB/sigmaS)*(s-my.muS)
    out<-pnorm(beta_0+beta_2*s+beta_4*c+beta_5)*pnorm( (c-mu3)/sigma3 ) * dnorm(s,mean=my.muS,sd=sigmaS) / pnorm( (c-my.muB)/sigmaB )
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
  area2<-Vectorize(area)
  
  term2<-area2(myid)
  
  risk<-term1+term2
  rr<-ifelse(dat$Y==1,risk,1-risk)
  ans<-ans+sum(log(rr))
  }
  
  #delta_S=1;delta_B=0,Z=1
  #if S>c 
  dat<-data[data$delta_S==1&data$delta_B==0&data$Z==1&data$S>c,] #65
  if(nrow(dat)>0){
    mu3<-ifelse(dat$W==1,muB[1]+rho*(sigmaB/sigmaS) *(dat$S-muS[1]),
                ifelse(dat$W==2, muB[2]+rho*(sigmaB/sigmaS) *(dat$S-muS[2]),
                       ifelse(dat$W==3, muB[3]+rho*(sigmaB/sigmaS) *(dat$S-muS[3]), muB[4]+rho*(sigmaB/sigmaS) *(dat$S-muS[4]))))
    beta_5<-ifelse(dat$W==1,0,ifelse(dat$W==2, beta_52, ifelse(dat$W==3, beta_53, beta_54)))
    term1<- pnorm( (c-mu3)/sigma3 ) * pnorm(beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*c+beta_5)

    myid<-1:nrow(dat)
    integrand<-function(b,id){
      if(dat$W[id]==1) {my.mu3<-muB[1]+rho*(sigmaB/sigmaS) *(dat$S[id]-muS[1]); beta_5<-0}
      if(dat$W[id]==2) {my.mu3<-muB[2]+rho*(sigmaB/sigmaS) *(dat$S[id]-muS[2]); beta_5<-beta_52}
      if(dat$W[id]==3) {my.mu3<-muB[3]+rho*(sigmaB/sigmaS) *(dat$S[id]-muS[3]); beta_5<-beta_53}
      if(dat$W[id]==4) {my.mu3<-muB[4]+rho*(sigmaB/sigmaS) *(dat$S[id]-muS[4]); beta_5<-beta_54}
      out<-pnorm(beta_0+beta_1+beta_2*dat$S[id]+beta_3*dat$S[id]+beta_4*b+beta_5) * dnorm(b,mean=my.mu3,sd=sigma3)
      return(out)
    }
    
    area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
    area2<-Vectorize(area)
    
    term2<-area2(myid)
    
    risk<-term1+term2
    rr<-ifelse(dat$Y==1,risk,1-risk)
    ans<-ans+sum(log(rr))
  }
  
  #if S=c
  dat<-data[data$delta_S==1&data$delta_B==0&data$Z==1&data$S==c,] #20
  if(nrow(dat)>0){
    p2<-as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[1],muS[1]),sigma = sigma))/pnorm( (c-muS[1])/sigmaS)
    p2<-c(p2,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[2],muS[2]),sigma = sigma))/pnorm( (c-muS[2])/sigmaS))
    p2<-c(p2,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[3],muS[3]),sigma = sigma))/pnorm( (c-muS[3])/sigmaS))
    p2<-c(p2,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[4],muS[4]),sigma = sigma))/pnorm( (c-muS[4])/sigmaS))
    
    term1<-ifelse(dat$W==1, p2[1]*pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*c),
                  ifelse(dat$W==2, p2[2]*pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*c+beta_52),
                         ifelse(dat$W==3,p2[3]*pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*c+beta_53),p2[4]*pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*c+beta_54))))
    
    myid<-1:nrow(dat)
    integrand<-function(b,id){
      if(dat$W[id]==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.muS<-muS[1]; my.muB<-muB[1]; beta_5<-0}
      if(dat$W[id]==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.muS<-muS[2]; my.muB<-muB[2]; beta_5<-beta_52}
      if(dat$W[id]==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.muS<-muS[3]; my.muB<-muB[3]; beta_5<-beta_53}
      if(dat$W[id]==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.muS<-muS[4]; my.muB<-muB[4]; beta_5<-beta_54}
      
      out<-pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*b+beta_5)*pnorm( (c-my.mu2)/sigma2 ) * dnorm(b,mean=my.muB,sd=sigmaB) / pnorm( (c-my.muS)/sigmaS)
      return(out)
    }
    
    area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
    area2<-Vectorize(area)
    
    term2<-area2(myid)
    
    risk<-term1+term2
    rr<-ifelse(dat$Y==1,risk,1-risk)
    ans<-ans+sum(log(rr))
  }
  
  #delta_S=1,delta_B=0,Z=0 or delta_S=0
  dat<-data[(data$delta_S==1&data$delta_B==0&data$Z==0)|data$delta_S==0,] #6383
  
  p3<-as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[1],muS[1]),sigma = sigma))
  p3<-c(p3,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[2],muS[2]),sigma = sigma)))
  p3<-c(p3,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[3],muS[3]),sigma = sigma)))
  p3<-c(p3,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[4],muS[4]),sigma = sigma)))
  
  term1<-ifelse(dat$W==1,p3[1]*pnorm(beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_4*c),
                ifelse(dat$W==2, p3[2]*pnorm(beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_4*c+beta_52),
                       ifelse(dat$W==3, p3[3]*pnorm(beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_4*c+beta_53), p3[4]*pnorm(beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_4*c+beta_54 ))))
  
  myid<-1:nrow(dat)
  integrand<-function(b,id){
    if(dat$W[id]==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]; beta_5<-0}
    if(dat$W[id]==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]; beta_5<-beta_52}
    if(dat$W[id]==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]; beta_5<-beta_53}
    if(dat$W[id]==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]; beta_5<-beta_54}
    
    out<-pnorm( (c-my.mu2)/sigma2 )* pnorm(beta_0+beta_1*dat$Z[id]+beta_2*c+beta_3*dat$Z[id]*c+beta_4*b+beta_5) * dnorm(b,mean=my.mu1,sd=sigma1)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
  area2<-Vectorize(area)
  term2<-area2(myid)
  
  integrand<-function(s,id){
    if(dat$W[id]==1) {my.mu3<-muB[1]+rho*(sigmaB/sigmaS) *(s-muS[1]); my.muS<-muS[1]; beta_5<-0}
    if(dat$W[id]==2) {my.mu3<-muB[2]+rho*(sigmaB/sigmaS) *(s-muS[2]); my.muS<-muS[2]; beta_5<-beta_52}
    if(dat$W[id]==3) {my.mu3<-muB[3]+rho*(sigmaB/sigmaS) *(s-muS[3]); my.muS<-muS[3]; beta_5<-beta_53}
    if(dat$W[id]==4) {my.mu3<-muB[4]+rho*(sigmaB/sigmaS) *(s-muS[4]); my.muS<-muS[4]; beta_5<-beta_54}
    
    out<-pnorm( (c-my.mu3)/sigma3 )* pnorm(beta_0+beta_1*dat$Z[id]+beta_2*s+beta_3*dat$Z[id]*s+beta_4*c+beta_5) * dnorm(s,mean=my.muS,sd=sigmaS)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
  area2<-Vectorize(area)
  term3<-area2(myid)
  
  #term4 only takes on 8 values, 4 (for 4 W values) for Z=1 and 4 for Z=0. 
  integrand<-function(s,b){
    mu2<-gamma0+gamma1*b+mygamma2
    out<-pnorm(beta_0+beta_1*1+beta_2*s+beta_3*1*s+beta_4*b+beta_5)* dnorm(s, mean=mu2,sd=sigma2) * dnorm(b,mean=my.mu1,sd=sigma1)
    return(out)
  }
  
  mygamma2<-gamma2[1]; my.mu1<-alpha0+alpha1[1]; beta_5<-0
  term4.1.1<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[2]; my.mu1<-alpha0+alpha1[2]; beta_5<-beta_52
  term4.1.2<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[3]; my.mu1<-alpha0+alpha1[3]; beta_5<-beta_53
  term4.1.3<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[4]; my.mu1<-alpha0+alpha1[4]; beta_5<-beta_54
  term4.1.4<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  
  integrand<-function(s,b){
    mu2<-gamma0+gamma1*b+mygamma2
    out<-pnorm(beta_0+beta_1*0+beta_2*s+beta_3*0*s+beta_4*b+beta_5)* dnorm(s, mean=mu2,sd=sigma2) * dnorm(b,mean=my.mu1,sd=sigma1)
    return(out)
  }
  
  mygamma2<-gamma2[1]; my.mu1<-alpha0+alpha1[1]; beta_5<-0
  term4.0.1<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[2]; my.mu1<-alpha0+alpha1[2]; beta_5<-beta_52
  term4.0.2<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[3]; my.mu1<-alpha0+alpha1[3]; beta_5<-beta_53
  term4.0.3<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[4]; my.mu1<-alpha0+alpha1[4]; beta_5<-beta_54
  term4.0.4<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  term4<-ifelse(dat$Z==1&dat$W==1,term4.1.1,
                ifelse(dat$Z==1&dat$W==2,term4.1.2,
                       ifelse(dat$Z==1&dat$W==3,term4.1.3,
                              ifelse(dat$Z==1&dat$W==4,term4.1.4,
                                     ifelse(dat$Z==0&dat$W==1,term4.0.1,
                                            ifelse(dat$Z==0&dat$W==2,term4.0.2,
                                                   ifelse(dat$Z==0&dat$W==3,term4.0.3,term4.0.4)))))))
  

 
  
#The Following lines is an example on how to do double integral for vector input/output, though it is slower than the for loop.  
   
#   myid<-1:nrow(dat)
#   integrand<-function(s,b,id){
#     if(dat$W[id]==1) {mygamma2<-gamma2[1]; my.mu1<-alpha0+alpha1[1]; beta_5<-0}
#     if(dat$W[id]==2) {mygamma2<-gamma2[2]; my.mu1<-alpha0+alpha1[2]; beta_5<-beta_52}
#     if(dat$W[id]==3) {mygamma2<-gamma2[3]; my.mu1<-alpha0+alpha1[3]; beta_5<-beta_53}
#     if(dat$W[id]==4) {mygamma2<-gamma2[4]; my.mu1<-alpha0+alpha1[4]; beta_5<-beta_54}
#     mu2<-gamma0+gamma1*b+mygamma2
#     out<-pnorm(beta_0+beta_1*dat$Z[id]+beta_2*s+beta_3*dat$Z[id]*s+beta_4*b+beta_5)* dnorm(s, mean=mu2,sd=sigma2) * dnorm(b,mean=my.mu1,sd=sigma1)
#     return(out)
#   }
#   
#   area<-function(id) integrate(function(b,id) { 
#     sapply(b, function(b) {
#       integrate(function(s){integrand(s,b,id)}, c, Inf)$value
#     })
#   }, c, Inf, id=id)$value
#     
#   area2<-Vectorize(area)
#   term4<-area2(myid)
    
    
  risk<-term1+term2+term3+term4
  rr<-ifelse(dat$Y==1,risk,1-risk)
  ans<-ans+sum(log(rr))
  
  #(Negative) log-likelihood function
  ans<- -ans
  
}




getBeta<-function(data, beta.start){
  ################Estimate Nuisance parameters
  #Estimate B|W
  # Estimate theta
  theta.start <- c(2,0.4,-0.6,0.3, 0.6)
  out1 <- optim(par=theta.start, fn=logl.BonW,data=data)$par
  #out<- nlm(logl.BonW, theta.start, data = data) produces the same results
  alpha0<-out1[1]
  alpha1<-c(0,out1[2:4])
  sigma1<-out1[5]

  #Estimate S|B,W
  #theta=c(gamma0,gamma1,gamma22,gamma23,gamma24,sigma2)
  # Estimate theta
  theta.start <- c(nuisance_true[[4]],nuisance_true[[5]],nuisance_true[[6]][-1], 0.6)
  out2 <- optim(par=theta.start, fn=logl.SonBW, data=data, alpha0=alpha0, alpha1=alpha1, sigma1=sigma1)$par
  gamma0<-out2[1]
  gamma1<-out2[2]
  gamma2<-c(0,out2[3:5])
  sigma2<-out2[6]
  
  muB<-alpha0+alpha1
  muS<-c(gamma0+gamma1*(alpha0+alpha1[1])+gamma2[1],
         gamma0+gamma1*(alpha0+alpha1[2])+gamma2[2],
         gamma0+gamma1*(alpha0+alpha1[3])+gamma2[3],
         gamma0+gamma1*(alpha0+alpha1[4])+gamma2[4])
  sigmaB<-sigma1
  sigmaS<- sqrt(sigma2^2+ (gamma1^2)*(sigma1^2))
  rho<- gamma1*sigma1/sqrt(sigma2^2+gamma1^2*sigma1^2)
  
  nuisance<-list(alpha0,alpha1, sigma1, gamma0, gamma1, gamma2, sigma2, muB, muS, sigmaB, sigmaS, rho)
  
  # Estimate beta
  out <- optim(par=beta.start, fn=like,data=data,nuisance=nuisance)$par
  out<-list(nuisance=nuisance, beta=out)
  return(out)
}



#################perturbed version

logl.BonW.Perturb <- function(theta,data,xi) {
  dat<-data[data$delta_B==1,]
  alpha0<-theta[1]
  alpha1<-c(0,theta[2:4])
  sigma1<-theta[5]
  mu1<-ifelse(dat$W==1,alpha0+alpha1[1],ifelse(dat$W==2,alpha0+alpha1[2],ifelse(dat$W==3,alpha0+alpha1[3],alpha0+alpha1[4])))
  ff<-ifelse(dat$B==c,pnorm(c,mean=mu1,sd=sigma1), dnorm(dat$B,mean=mu1,sd=sigma1))
  ff<-log(ff)*xi[data$delta_B==1]
  return(-sum(ff))
}


logl.SonBW.Perturb<-function(theta,data,alpha0,alpha1,sigma1,xi){
  dat<-data[data$delta_B==1&data$Z==1,]
  dat$xi<-xi[data$delta_B==1&data$Z==1]
  gamma0<-theta[1]
  gamma1<-theta[2]
  gamma2<-c(0,theta[3:5])
  sigma2<-theta[6]
  
  ans<-0
  #S=c,B=c
  dd<-dat[dat$S==c&dat$B==c,] #6
  if(nrow(dd)>0){
  myid<-1:nrow(dd)
  integrand<-function(b,id){
    if(dd$W[id]==1) {mygamma2<-gamma2[1]; mu1<-alpha0+alpha1[1]}
    if(dd$W[id]==2) {mygamma2<-gamma2[2]; mu1<-alpha0+alpha1[2]}
    if(dd$W[id]==3) {mygamma2<-gamma2[3]; mu1<-alpha0+alpha1[3]}
    if(dd$W[id]==4) {mygamma2<-gamma2[4]; mu1<-alpha0+alpha1[4]}
    
    out<-pnorm(c, mean= gamma0+gamma1*b+mygamma2,sd=sigma2)*dnorm(b, mean=mu1,sd=sigma1)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=-Inf,upper=c,id=id)$value
  area2<-Vectorize(area)
  num<-area2(myid)
  
  mu1<-ifelse(dd$W==1,alpha0+alpha1[1],ifelse(dd$W==2,alpha0+alpha1[2],ifelse(dd$W==3,alpha0+alpha1[3],alpha0+alpha1[4])))
  den<-pnorm(c,mean=mu1,sd=sigma1)
  ff<-num/den
  ff<-log(ff)*dd$xi
  ans<-ans+sum(ff)
  }
  
  #S=c,B>c
  dd<-dat[dat$S==c&dat$B>c,] #0
  if(nrow(dd)>0){
    mu2<-ifelse(dd$W==1,gamma0+gamma1*dd$B+gamma2[1],ifelse(dd$W==2,gamma0+gamma1*dd$B+gamma2[2],ifelse(dd$W==3,gamma0+gamma1*dd$B+gamma2[3],gamma0+gamma1*dd$B+gamma2[4])))
    ff<-dnorm(c,mean=mu2,sd=sigma2)
    ff<-log(ff)*dd$xi
    ans<-ans+sum(ff)
  }
  
  #S>c,B=c
  dd<-dat[dat$S>c&dat$B==c,] #101
  if(nrow(dd)>0){
  myid<-1:nrow(dd)
  integrand<-function(b,id){
    if(dd$W[id]==1) {mygamma2<-gamma2[1]; mu1<-alpha0+alpha1[1]}
    if(dd$W[id]==2) {mygamma2<-gamma2[2]; mu1<-alpha0+alpha1[2]}
    if(dd$W[id]==3) {mygamma2<-gamma2[3]; mu1<-alpha0+alpha1[3]}
    if(dd$W[id]==4) {mygamma2<-gamma2[4]; mu1<-alpha0+alpha1[4]}
    
    out<-dnorm(dd$S[id], mean= gamma0+gamma1*b+mygamma2,sd=sigma2)*dnorm(b, mean=mu1,sd=sigma1)
    return(out)
  }
  area<-function(id) integrate(integrand,lower=-Inf,upper=c,id=id)$value
  area2<-Vectorize(area)
  num<-area2(myid)
  
  mu1<-ifelse(dd$W==1,alpha0+alpha1[1],ifelse(dd$W==2,alpha0+alpha1[2],ifelse(dd$W==3,alpha0+alpha1[3],alpha0+alpha1[4])))
  den<-pnorm(c,mean=mu1,sd=sigma1)
  ff<-num/den
  ff<-log(ff)*dd$xi
  ans<-ans+sum(ff)
  }
  
  #S>c,B>c
  dd<-dat[dat$S>c&dat$B>c,] #2304
  mu2<-ifelse(dd$W==1,gamma0+gamma1*dd$B+gamma2[1],ifelse(dd$W==2,gamma0+gamma1*dd$B+gamma2[2],ifelse(dd$W==3,gamma0+gamma1*dd$B+gamma2[3],gamma0+gamma1*dd$B+gamma2[4])))
  ff<-dnorm(dd$S,mean=mu2,sd=sigma2)
  ff<-log(ff)*dd$xi
  ans<-ans+sum(ff)
  
  return(-ans)
}





like.Perturb<-function(beta,data,nuisance,xi){
  beta_0<-beta[1];beta_1<-beta[2];beta_2<-beta[3];beta_3<-beta[4];beta_4<-beta[5];beta_52<-beta[6];beta_53<-beta[7];beta_54<-beta[8]
  
  alpha0<-nuisance[[1]]
  alpha1<-nuisance[[2]]
  sigma1<-nuisance[[3]]
  gamma0<-nuisance[[4]]
  gamma1<-nuisance[[5]]
  gamma2<-nuisance[[6]]
  sigma2<-nuisance[[7]]
  muB<-nuisance[[8]]
  muS<-nuisance[[9]]
  sigmaB<-nuisance[[10]]
  sigmaS<-nuisance[[11]]
  rho<-nuisance[[12]]
  sigma3<-sigmaB*sqrt(1-rho^2)
  sigma<-matrix(c(sigmaB^2,rho*sigmaB*sigmaS,rho*sigmaB*sigmaS,sigmaS^2),nrow=2)
  
  ans <- 0
  data$xi<-xi
  # delta_S=delta_B=1,Z=1
  dat<-data[data$delta_S==1&data$delta_B==1&data$Z==1,] #2415
  ff<-ifelse(dat$W==1,beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*dat$B,
             ifelse(dat$W==2,beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*dat$B+beta_52,
                    ifelse(dat$W==3, beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*dat$B+beta_53, beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*dat$B+beta_54)))
  
  risk<-pnorm(ff)
  rr<-ifelse(dat$Y==1,risk,1-risk)
  rr<-log(rr)*dat$xi
  ans<-ans+sum(rr)
  
  #delta_S=delta_B=1,Z=0
  #if B>c 
  dat<-data[data$delta_S==1&data$delta_B==1&data$Z==0&data$B>c,] #1053
  mu2<- ifelse(dat$W==1, gamma0+gamma1*dat$B+gamma2[1],
               ifelse(dat$W==2, gamma0+gamma1*dat$B+gamma2[2],
                      ifelse(dat$W==3, gamma0+gamma1*dat$B+gamma2[3],gamma0+gamma1*dat$B+gamma2[4])))
  ff<-ifelse(dat$W==1,beta_0+beta_2*c+beta_4*dat$B,
             ifelse(dat$W==2,beta_0+beta_2*c+beta_4*dat$B+beta_52,
                    ifelse(dat$W==3, beta_0+beta_2*c+beta_4*dat$B+beta_53, beta_0+beta_2*c+beta_4*dat$B+beta_54)))
  term1<- pnorm( (c-mu2)/sigma2 ) * pnorm(ff)
  
  myid<-1:nrow(dat)
  integrand<-function(s,id){
    if(dat$W[id]==1) {mu2<-gamma0+gamma1*dat$B[id]+gamma2[1]; beta5<-0}
    if(dat$W[id]==2) {mu2<-gamma0+gamma1*dat$B[id]+gamma2[2]; beta5<-beta_52}
    if(dat$W[id]==3) {mu2<-gamma0+gamma1*dat$B[id]+gamma2[3]; beta5<-beta_53}
    if(dat$W[id]==4) {mu2<-gamma0+gamma1*dat$B[id]+gamma2[4]; beta5<-beta_54}
    
    out<-pnorm(beta_0+beta_2*s+beta_4*dat$B[id]+beta5)*dnorm(s, mean=mu2, sd=sigma2)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
  area2<-Vectorize(area)
  
  term2<-area2(myid)
  
  risk<-term1+term2
  rr<-ifelse(dat$Y==1,risk,1-risk)
  rr<-log(rr)*dat$xi
  ans<-ans+sum(rr)
  
  #if B=c
  dat<-data[data$delta_S==1&data$delta_B==1&data$Z==0&data$B==c,] #64
  if(nrow(dat)>0){
    p1<-as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[1],muS[1]),sigma = sigma))/pnorm( (c-muB[1])/sigmaB)
    p1<-c(p1,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[2],muS[2]),sigma = sigma))/pnorm( (c-muB[2])/sigmaB))
    p1<-c(p1,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[3],muS[3]),sigma = sigma))/pnorm( (c-muB[3])/sigmaB))
    p1<-c(p1,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[4],muS[4]),sigma = sigma))/pnorm( (c-muB[4])/sigmaB))
    
    term1<-ifelse(dat$W==1, p1[1]*pnorm(beta_0+beta_2*c+beta_4*c),
                  ifelse(dat$W==2, p1[2]*pnorm(beta_0+beta_2*c+beta_4*c+beta_52),
                         ifelse(dat$W==3,p1[3]*pnorm(beta_0+beta_2*c+beta_4*c+beta_53),p1[4]*pnorm(beta_0+beta_2*c+beta_4*c+beta_54))))
    
    myid<-1:nrow(dat)
    integrand<-function(s,id){
      if(dat$W[id]==1) {my.muS<-muS[1]; my.muB<-muB[1]; beta_5<-0}
      if(dat$W[id]==2) {my.muS<-muS[2]; my.muB<-muB[2]; beta_5<-beta_52}
      if(dat$W[id]==3) {my.muS<-muS[3]; my.muB<-muB[3]; beta_5<-beta_53}
      if(dat$W[id]==4) {my.muS<-muS[4]; my.muB<-muB[4]; beta_5<-beta_54}
      
      mu3<-my.muB+rho*(sigmaB/sigmaS)*(s-my.muS)
      out<-pnorm(beta_0+beta_2*s+beta_4*c+beta_5)*pnorm( (c-mu3)/sigma3 ) * dnorm(s,mean=my.muS,sd=sigmaS) / pnorm( (c-my.muB)/sigmaB )
      return(out)
    }
    
    area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
    area2<-Vectorize(area)
    
    term2<-area2(myid)
    
    risk<-term1+term2
    rr<-ifelse(dat$Y==1,risk,1-risk)
    rr<-log(rr)*dat$xi
    ans<-ans+sum(rr)
  }
  
  #delta_S=1;delta_B=0,Z=1
  #if S>c 
  dat<-data[data$delta_S==1&data$delta_B==0&data$Z==1&data$S>c,] #65
  if(nrow(dat)>0){
    mu3<-ifelse(dat$W==1,muB[1]+rho*(sigmaB/sigmaS) *(dat$S-muS[1]),
                ifelse(dat$W==2, muB[2]+rho*(sigmaB/sigmaS) *(dat$S-muS[2]),
                       ifelse(dat$W==3, muB[3]+rho*(sigmaB/sigmaS) *(dat$S-muS[3]), muB[4]+rho*(sigmaB/sigmaS) *(dat$S-muS[4]))))
    beta_5<-ifelse(dat$W==1,0,ifelse(dat$W==2, beta_52, ifelse(dat$W==3, beta_53, beta_54)))
    term1<- pnorm( (c-mu3)/sigma3 ) * pnorm(beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_4*c+beta_5)
    
    myid<-1:nrow(dat)
    integrand<-function(b,id){
      if(dat$W[id]==1) {my.mu3<-muB[1]+rho*(sigmaB/sigmaS) *(dat$S[id]-muS[1]); beta_5<-0}
      if(dat$W[id]==2) {my.mu3<-muB[2]+rho*(sigmaB/sigmaS) *(dat$S[id]-muS[2]); beta_5<-beta_52}
      if(dat$W[id]==3) {my.mu3<-muB[3]+rho*(sigmaB/sigmaS) *(dat$S[id]-muS[3]); beta_5<-beta_53}
      if(dat$W[id]==4) {my.mu3<-muB[4]+rho*(sigmaB/sigmaS) *(dat$S[id]-muS[4]); beta_5<-beta_54}
      out<-pnorm(beta_0+beta_1+beta_2*dat$S[id]+beta_3*dat$S[id]+beta_4*b+beta_5) * dnorm(b,mean=my.mu3,sd=sigma3)
      return(out)
    }
    
    area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
    area2<-Vectorize(area)
    
    term2<-area2(myid)
    
    risk<-term1+term2
    rr<-ifelse(dat$Y==1,risk,1-risk)
    rr<-log(rr)*dat$xi
    ans<-ans+sum(rr)
  }
  
  #if S=c
  dat<-data[data$delta_S==1&data$delta_B==0&data$Z==1&data$S==c,] #20
  if(nrow(dat)>0){
    p2<-as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[1],muS[1]),sigma = sigma))/pnorm( (c-muS[1])/sigmaS)
    p2<-c(p2,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[2],muS[2]),sigma = sigma))/pnorm( (c-muS[2])/sigmaS))
    p2<-c(p2,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[3],muS[3]),sigma = sigma))/pnorm( (c-muS[3])/sigmaS))
    p2<-c(p2,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[4],muS[4]),sigma = sigma))/pnorm( (c-muS[4])/sigmaS))
    
    term1<-ifelse(dat$W==1, p2[1]*pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*c),
                  ifelse(dat$W==2, p2[2]*pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*c+beta_52),
                         ifelse(dat$W==3,p2[3]*pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*c+beta_53),p2[4]*pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*c+beta_54))))
    
    myid<-1:nrow(dat)
    integrand<-function(b,id){
      if(dat$W[id]==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.muS<-muS[1]; my.muB<-muB[1]; beta_5<-0}
      if(dat$W[id]==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.muS<-muS[2]; my.muB<-muB[2]; beta_5<-beta_52}
      if(dat$W[id]==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.muS<-muS[3]; my.muB<-muB[3]; beta_5<-beta_53}
      if(dat$W[id]==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.muS<-muS[4]; my.muB<-muB[4]; beta_5<-beta_54}
      
      out<-pnorm(beta_0+beta_1+beta_2*c+beta_3*c+beta_4*b+beta_5)*pnorm( (c-my.mu2)/sigma2 ) * dnorm(b,mean=my.muB,sd=sigmaB) / pnorm( (c-my.muS)/sigmaS)
      return(out)
    }
    
    area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
    area2<-Vectorize(area)
    
    term2<-area2(myid)
    
    risk<-term1+term2
    rr<-ifelse(dat$Y==1,risk,1-risk)
    rr<-log(rr)*dat$xi
    ans<-ans+sum(rr)
  }
  
  #delta_S=1,delta_B=0,Z=0 or delta_S=0
  dat<-data[(data$delta_S==1&data$delta_B==0&data$Z==0)|data$delta_S==0,] #6383
  
  p3<-as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[1],muS[1]),sigma = sigma))
  p3<-c(p3,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[2],muS[2]),sigma = sigma)))
  p3<-c(p3,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[3],muS[3]),sigma = sigma)))
  p3<-c(p3,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[4],muS[4]),sigma = sigma)))
  
  term1<-ifelse(dat$W==1,p3[1]*pnorm(beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_4*c),
                ifelse(dat$W==2, p3[2]*pnorm(beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_4*c+beta_52),
                       ifelse(dat$W==3, p3[3]*pnorm(beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_4*c+beta_53), p3[4]*pnorm(beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_4*c+beta_54 ))))
  
  myid<-1:nrow(dat)
  integrand<-function(b,id){
    if(dat$W[id]==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]; beta_5<-0}
    if(dat$W[id]==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]; beta_5<-beta_52}
    if(dat$W[id]==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]; beta_5<-beta_53}
    if(dat$W[id]==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]; beta_5<-beta_54}
    
    out<-pnorm( (c-my.mu2)/sigma2 )* pnorm(beta_0+beta_1*dat$Z[id]+beta_2*c+beta_3*dat$Z[id]*c+beta_4*b+beta_5) * dnorm(b,mean=my.mu1,sd=sigma1)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
  area2<-Vectorize(area)
  term2<-area2(myid)
  
  integrand<-function(s,id){
    if(dat$W[id]==1) {my.mu3<-muB[1]+rho*(sigmaB/sigmaS) *(s-muS[1]); my.muS<-muS[1]; beta_5<-0}
    if(dat$W[id]==2) {my.mu3<-muB[2]+rho*(sigmaB/sigmaS) *(s-muS[2]); my.muS<-muS[2]; beta_5<-beta_52}
    if(dat$W[id]==3) {my.mu3<-muB[3]+rho*(sigmaB/sigmaS) *(s-muS[3]); my.muS<-muS[3]; beta_5<-beta_53}
    if(dat$W[id]==4) {my.mu3<-muB[4]+rho*(sigmaB/sigmaS) *(s-muS[4]); my.muS<-muS[4]; beta_5<-beta_54}
    
    out<-pnorm( (c-my.mu3)/sigma3 )* pnorm(beta_0+beta_1*dat$Z[id]+beta_2*s+beta_3*dat$Z[id]*s+beta_4*c+beta_5) * dnorm(s,mean=my.muS,sd=sigmaS)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
  area2<-Vectorize(area)
  term3<-area2(myid)
  
  #term4 only takes on 8 values, 4 (for 4 W values) for Z=1 and 4 for Z=0. 
  integrand<-function(s,b){
    mu2<-gamma0+gamma1*b+mygamma2
    out<-pnorm(beta_0+beta_1*1+beta_2*s+beta_3*1*s+beta_4*b+beta_5)* dnorm(s, mean=mu2,sd=sigma2) * dnorm(b,mean=my.mu1,sd=sigma1)
    return(out)
  }
  
  mygamma2<-gamma2[1]; my.mu1<-alpha0+alpha1[1]; beta_5<-0
  term4.1.1<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[2]; my.mu1<-alpha0+alpha1[2]; beta_5<-beta_52
  term4.1.2<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[3]; my.mu1<-alpha0+alpha1[3]; beta_5<-beta_53
  term4.1.3<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[4]; my.mu1<-alpha0+alpha1[4]; beta_5<-beta_54
  term4.1.4<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  
  integrand<-function(s,b){
    mu2<-gamma0+gamma1*b+mygamma2
    out<-pnorm(beta_0+beta_1*0+beta_2*s+beta_3*0*s+beta_4*b+beta_5)* dnorm(s, mean=mu2,sd=sigma2) * dnorm(b,mean=my.mu1,sd=sigma1)
    return(out)
  }
  
  mygamma2<-gamma2[1]; my.mu1<-alpha0+alpha1[1]; beta_5<-0
  term4.0.1<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[2]; my.mu1<-alpha0+alpha1[2]; beta_5<-beta_52
  term4.0.2<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[3]; my.mu1<-alpha0+alpha1[3]; beta_5<-beta_53
  term4.0.3<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  mygamma2<-gamma2[4]; my.mu1<-alpha0+alpha1[4]; beta_5<-beta_54
  term4.0.4<-integrate(function(b) { 
    sapply(b, function(b) {
      integrate(function(s){integrand(s,b)}, c, Inf)$value
    })
  }, c, Inf)$value
  
  term4<-ifelse(dat$Z==1&dat$W==1,term4.1.1,
                ifelse(dat$Z==1&dat$W==2,term4.1.2,
                       ifelse(dat$Z==1&dat$W==3,term4.1.3,
                              ifelse(dat$Z==1&dat$W==4,term4.1.4,
                                     ifelse(dat$Z==0&dat$W==1,term4.0.1,
                                            ifelse(dat$Z==0&dat$W==2,term4.0.2,
                                                   ifelse(dat$Z==0&dat$W==3,term4.0.3,term4.0.4)))))))
  
  
  
  
  #The Following lines is an example on how to do double integral for vector input/output, though it is slower than the for loop.  
  
  #   myid<-1:nrow(dat)
  #   integrand<-function(s,b,id){
  #     if(dat$W[id]==1) {mygamma2<-gamma2[1]; my.mu1<-alpha0+alpha1[1]; beta_5<-0}
  #     if(dat$W[id]==2) {mygamma2<-gamma2[2]; my.mu1<-alpha0+alpha1[2]; beta_5<-beta_52}
  #     if(dat$W[id]==3) {mygamma2<-gamma2[3]; my.mu1<-alpha0+alpha1[3]; beta_5<-beta_53}
  #     if(dat$W[id]==4) {mygamma2<-gamma2[4]; my.mu1<-alpha0+alpha1[4]; beta_5<-beta_54}
  #     mu2<-gamma0+gamma1*b+mygamma2
  #     out<-pnorm(beta_0+beta_1*dat$Z[id]+beta_2*s+beta_3*dat$Z[id]*s+beta_4*b+beta_5)* dnorm(s, mean=mu2,sd=sigma2) * dnorm(b,mean=my.mu1,sd=sigma1)
  #     return(out)
  #   }
  #   
  #   area<-function(id) integrate(function(b,id) { 
  #     sapply(b, function(b) {
  #       integrate(function(s){integrand(s,b,id)}, c, Inf)$value
  #     })
  #   }, c, Inf, id=id)$value
  #     
  #   area2<-Vectorize(area)
  #   term4<-area2(myid)
  
  
  risk<-term1+term2+term3+term4
  rr<-ifelse(dat$Y==1,risk,1-risk)
  rr<-log(rr)*dat$xi
  ans<-ans+sum(rr)
  
  #(Negative) log-likelihood function
  ans<- -ans
  
}



getBeta.Perturb<-function(data, beta.start, xi){
  ################Estimate Nuisance parameters
  #Estimate B|W
  # Estimate theta
  theta.start <- c(2,0.4,-0.6,0.3, 0.6)
  out1 <- optim(par=theta.start, fn=logl.BonW.Perturb,data=data,xi=xi)$par
  #out<- nlm(logl.BonW, theta.start, data = data) produces the same results
  alpha0<-out1[1]
  alpha1<-c(0,out1[2:4])
  sigma1<-out1[5]
  
  #Estimate S|B,W
  #theta=c(gamma0,gamma1,gamma22,gamma23,gamma24,sigma2)
  # Estimate theta
  theta.start <- c(nuisance_true[[4]],nuisance_true[[5]],nuisance_true[[6]][-1], 0.6)
  out2 <- optim(par=theta.start, fn=logl.SonBW.Perturb, data=data, alpha0=alpha0, alpha1=alpha1, sigma1=sigma1, xi=xi)$par
  gamma0<-out2[1]
  gamma1<-out2[2]
  gamma2<-c(0,out2[3:5])
  sigma2<-out2[6]
  
  muB<-alpha0+alpha1
  muS<-c(gamma0+gamma1*(alpha0+alpha1[1])+gamma2[1],
         gamma0+gamma1*(alpha0+alpha1[2])+gamma2[2],
         gamma0+gamma1*(alpha0+alpha1[3])+gamma2[3],
         gamma0+gamma1*(alpha0+alpha1[4])+gamma2[4])
  sigmaB<-sigma1
  sigmaS<- sqrt(sigma2^2+ (gamma1^2)*(sigma1^2))
  rho<- gamma1*sigma1/sqrt(sigma2^2+gamma1^2*sigma1^2)
  
  nuisance<-list(alpha0,alpha1, sigma1, gamma0, gamma1, gamma2, sigma2, muB, muS, sigmaB, sigmaS, rho)
  
  # Estimate beta
  out <- optim(par=beta.start, fn=like.Perturb,data=data,nuisance=nuisance, xi=xi)$par
  out<-list(nuisance=nuisance, beta=out)
  return(out)
}

#risk_z(S(1),B,W) can be calculate using beta

#risk_z(S(1),W)
riskonSW<-function(s,w,z,nuisance,beta){
  beta_0<-beta[1];beta_1<-beta[2];beta_2<-beta[3];beta_3<-beta[4];beta_4<-beta[5];beta_52<-beta[6];beta_53<-beta[7];beta_54<-beta[8]
  
  alpha0<-nuisance[[1]]
  alpha1<-nuisance[[2]]
  sigma1<-nuisance[[3]]
  gamma0<-nuisance[[4]]
  gamma1<-nuisance[[5]]
  gamma2<-nuisance[[6]]
  sigma2<-nuisance[[7]]
  muB<-nuisance[[8]]
  muS<-nuisance[[9]]
  sigmaB<-nuisance[[10]]
  sigmaS<-nuisance[[11]]
  rho<-nuisance[[12]]
  sigma3<-sigmaB*sqrt(1-rho^2)
  sigma<-matrix(c(sigmaB^2,rho*sigmaB*sigmaS,rho*sigmaB*sigmaS,sigmaS^2),nrow=2)
  
  if(s>c){
    mu3<-ifelse(w==1,muB[1]+rho*(sigmaB/sigmaS) *(s-muS[1]),
                ifelse(w==2, muB[2]+rho*(sigmaB/sigmaS) *(s-muS[2]),
                       ifelse(w==3, muB[3]+rho*(sigmaB/sigmaS) *(s-muS[3]), muB[4]+rho*(sigmaB/sigmaS) *(s-muS[4]))))
    
    beta_5<-ifelse(w==1,0,ifelse(w==2, beta_52, ifelse(w==3, beta_53, beta_54)))
    term1<- pnorm( (c-mu3)/sigma3 ) * pnorm(beta_0+beta_1*z+beta_2*s+beta_3*z*s+beta_4*c+beta_5)
    
    integrand<-function(b){
      out<-pnorm(beta_0+beta_1*z+beta_2*s+beta_3*z*s+beta_4*b+beta_5) * dnorm(b,mean=mu3,sd=sigma3)
      return(out)
    }
    
    term2<-integrate(integrand,lower=c,upper=Inf)$value
    
    out<-term1+term2
  }
  
  if(s==c){
    p2<-ifelse(w==1,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[1],muS[1]),sigma = sigma))/pnorm( (c-muS[1])/sigmaS),
               ifelse(w==2,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[2],muS[2]),sigma = sigma))/pnorm( (c-muS[2])/sigmaS),
                      ifelse(w==3,as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[3],muS[3]),sigma = sigma))/pnorm( (c-muS[3])/sigmaS),
                             as.numeric(pmvnorm(lower=c(-Inf,-Inf), upper=c(c,c), mean=c(muB[4],muS[4]),sigma = sigma))/pnorm( (c-muS[4])/sigmaS))))
    
    beta_5<-ifelse(w==1,0,ifelse(w==2, beta_52, ifelse(w==3, beta_53, beta_54)))
    term1<-p2*pnorm(beta_0+beta_1*z+beta_2*c+beta_3*z*c+beta_4*c+beta_5)
    
    integrand<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.muS<-muS[1]; my.muB<-muB[1]; beta_5<-0}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.muS<-muS[2]; my.muB<-muB[2]; beta_5<-beta_52}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.muS<-muS[3]; my.muB<-muB[3]; beta_5<-beta_53}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.muS<-muS[4]; my.muB<-muB[4]; beta_5<-beta_54}
      
      out<-pnorm(beta_0+beta_1*z+beta_2*c+beta_3*z*c+beta_4*b+beta_5)*pnorm( (c-my.mu2)/sigma2 ) * dnorm(b,mean=my.muB,sd=sigmaB) / pnorm( (c-my.muS)/sigmaS)
      return(out)
    }
    term2<-integrate(integrand,lower=c,upper=Inf)$value
    
    out<-term1+term2
  }
  return(out)
}


#risk_z(S(1),W) on subset defined on B>c
riskonSW_Bpos<-function(s,w,z,nuisance,beta){
  beta_0<-beta[1];beta_1<-beta[2];beta_2<-beta[3];beta_3<-beta[4];beta_4<-beta[5];beta_52<-beta[6];beta_53<-beta[7];beta_54<-beta[8]
  
  alpha0<-nuisance[[1]]
  alpha1<-nuisance[[2]]
  sigma1<-nuisance[[3]]
  gamma0<-nuisance[[4]]
  gamma1<-nuisance[[5]]
  gamma2<-nuisance[[6]]
  sigma2<-nuisance[[7]]
  muB<-nuisance[[8]]
  muS<-nuisance[[9]]
  sigmaB<-nuisance[[10]]
  sigmaS<-nuisance[[11]]
  rho<-nuisance[[12]]
  sigma3<-sigmaB*sqrt(1-rho^2)
  sigma<-matrix(c(sigmaB^2,rho*sigmaB*sigmaS,rho*sigmaB*sigmaS,sigmaS^2),nrow=2)
  
  if(s>c){
    mu3<-ifelse(w==1,muB[1]+rho*(sigmaB/sigmaS) *(s-muS[1]),
                ifelse(w==2, muB[2]+rho*(sigmaB/sigmaS) *(s-muS[2]),
                       ifelse(w==3, muB[3]+rho*(sigmaB/sigmaS) *(s-muS[3]), muB[4]+rho*(sigmaB/sigmaS) *(s-muS[4]))))
    
    beta_5<-ifelse(w==1,0,ifelse(w==2, beta_52, ifelse(w==3, beta_53, beta_54)))
    
    integrand<-function(b){
      out<-pnorm(beta_0+beta_1*z+beta_2*s+beta_3*z*s+beta_4*b+beta_5) * dnorm(b,mean=mu3,sd=sigma3)
      return(out)
    }
    
    num<-integrate(integrand,lower=c,upper=Inf)$value
    den<-1-pnorm(c,mean=mu3,sd=sigma3)
    out<-num/den
  }
  
  if(s==c){
    
    ff<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.muS<-muS[1]; my.muB<-muB[1]; beta_5<-0}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.muS<-muS[2]; my.muB<-muB[2]; beta_5<-beta_52}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.muS<-muS[3]; my.muB<-muB[3]; beta_5<-beta_53}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.muS<-muS[4]; my.muB<-muB[4]; beta_5<-beta_54}
      
      out<-pnorm( (c-my.mu2)/sigma2 ) * dnorm(b,mean=my.muB,sd=sigmaB) / pnorm( (c-my.muS)/sigmaS)
      return(out)
    }
    den<-integrate(ff,lower=c,upper=Inf)$value
    
    integrand<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.muS<-muS[1]; my.muB<-muB[1]; beta_5<-0}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.muS<-muS[2]; my.muB<-muB[2]; beta_5<-beta_52}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.muS<-muS[3]; my.muB<-muB[3]; beta_5<-beta_53}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.muS<-muS[4]; my.muB<-muB[4]; beta_5<-beta_54}
      
      out<-pnorm(beta_0+beta_1*z+beta_2*s+beta_3*z*s+beta_4*b+beta_5)*ff(b)
      return(out)
    }
    
    num<-integrate(integrand,lower=c,upper=Inf)$value
    out<-num/den
  }
  return(out)
}


#risk_z(S(1),W) on subset defined on B=c
riskonSW_Bneg<-function(s,w,z,nuisance,beta){
  beta_0<-beta[1];beta_1<-beta[2];beta_2<-beta[3];beta_3<-beta[4];beta_4<-beta[5];beta_52<-beta[6];beta_53<-beta[7];beta_54<-beta[8]
  if(w==1) {beta_5<-0}
  if(w==2) {beta_5<-beta_52}
  if(w==3) {beta_5<-beta_53}
  if(w==4) {beta_5<-beta_54}
  out<-pnorm(beta_0+beta_1*z+beta_2*s+beta_3*z*s+beta_4*c+beta_5)
  return(out)
}




#######get VE:
getVE<-function(data,Su,results,type){
  nuisance<-results$nuisance
  alpha0<-nuisance[[1]]
  alpha1<-nuisance[[2]]
  sigma1<-nuisance[[3]]
  gamma0<-nuisance[[4]]
  gamma1<-nuisance[[5]]
  gamma2<-nuisance[[6]]
  sigma2<-nuisance[[7]]
  muB<-nuisance[[8]]
  muS<-nuisance[[9]]
  sigmaB<-nuisance[[10]]
  sigmaS<-nuisance[[11]]
  rho<-nuisance[[12]]
  sigma3<-sigmaB*sqrt(1-rho^2)
  sigma<-matrix(c(sigmaB^2,rho*sigmaB*sigmaS,rho*sigmaB*sigmaS,sigmaS^2),nrow=2)
   
  data$S1<-ifelse(data$Z==1,data$S,NA)
  Wu<-sort(unique(data$W[!is.na(data$W)]))
  
  if(type=="marginal"){
    myriskonSW<-riskonSW
    fX<-table(data$W)/sum(!is.na(data$W))
    #When s>c
    integrand<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- dnorm(mys,mean=my.mu2,sd=sigma2)*dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    
    fx.s<-NULL
    for(mys in Su[-1]){
      myrow<-NULL
      for(w in 1:4){
    myrow<-c(myrow,integrate(integrand,lower=-Inf,upper=Inf)$value  )
      }
    myrow<-myrow*fX/sum(myrow*fX)
    fx.s<-rbind(fx.s,myrow)
    }
    
    #When s=c
    integrand<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- pnorm(c,mean=my.mu2,sd=sigma2)*dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    myrow<-NULL
    for(w in 1:4){
      myrow<-c(myrow,integrate(integrand,lower=-Inf,upper=Inf)$value  )
    }
    myrow<-myrow*fX/sum(myrow*fX)
    fx.s<-rbind(myrow,fx.s)
  }
  
  if(type=="Bpos"){
    myriskonSW<-riskonSW_Bpos
    #P(W|B>c)
    fX<-table(data$W)/sum(!is.na(data$W))
    myrow<-NULL
    for(w in 1:4){
      if(w==1) {my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu1<-alpha0+alpha1[4]}
      myrow<-c(myrow,1-pnorm(c,mean=my.mu1,sd=sigma1))
    }
    fX<-myrow*fX/sum(myrow*fX)
    
    #When s>c
    integrand1<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- dnorm(mys,mean=my.mu2,sd=sigma2)*dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    integrand2<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
  
    fx.s<-NULL
    for(mys in Su[-1]){
      myrow<-NULL
      for(w in 1:4){
        myrow<-c(myrow,integrate(integrand1,lower=c,upper=Inf)$value/integrate(integrand2,lower=c,upper=Inf)$value  )
      }
      myrow<-myrow*fX/sum(myrow*fX)
      fx.s<-rbind(fx.s,myrow)
    }
    
    #When s=c
    integrand1<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- pnorm(c,mean=my.mu2,sd=sigma2)*dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    myrow<-NULL
    for(w in 1:4){
      myrow<-c(myrow,integrate(integrand1,lower=c,upper=Inf)$value/integrate(integrand2,lower=c,upper=Inf)$value  )
    }
    myrow<-myrow*fX/sum(myrow*fX)
    fx.s<-rbind(myrow,fx.s)
  }
  
  if(type=="Bneg"){
    myriskonSW<-riskonSW_Bneg
    #P(W|B=c)
    fX<-table(data$W)/sum(!is.na(data$W))
    myrow<-NULL
    for(w in 1:4){
      if(w==1) {my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu1<-alpha0+alpha1[4]}
      myrow<-c(myrow,pnorm(c,mean=my.mu1,sd=sigma1))
    }
    fX<-myrow*fX/sum(myrow*fX)
    
    #When s>c
    b<-c
    fx.s<-NULL
    for(mys in Su[-1]){
      myrow<-NULL
      for(w in 1:4){
        if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
        if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
        if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
        if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
        myrow<-c(myrow, dnorm(mys,mean=my.mu2,sd=sigma2) )
      }
      myrow<-myrow*fX/sum(myrow*fX)
      fx.s<-rbind(fx.s,myrow)
    }
    #When s=c
    b<-c
    mys <-c
    myrow<-NULL
      for(w in 1:4){
        if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
        if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
        if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
        if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
        myrow<-c(myrow, pnorm(mys,mean=my.mu2,sd=sigma2) )
      }
      myrow<-myrow*fX/sum(myrow*fX)
      fx.s<-rbind(myrow,fx.s)
    }
  
  PY1.SZ0<-PY1.SZ1<-rep(0,length(Su))
  
  for(i in 1:length(Su)){
    for (j in 1:length(Wu)){
      PY1.SZ1[i]=PY1.SZ1[i]+myriskonSW(s=Su[i],w=Wu[j],z=1,nuisance = results$nuisance, beta = results$beta)*fx.s[i,j]
      PY1.SZ0[i]=PY1.SZ0[i]+myriskonSW(s=Su[i],w=Wu[j],z=0,nuisance = results$nuisance, beta = results$beta)*fx.s[i,j]
    }
  }
  #in case PY1.SZ0 is very very small, and R recorded as 0:
  PY1.SZ0<-ifelse(PY1.SZ0==0,1e-300,PY1.SZ0)
  VE.S<-1-PY1.SZ1/PY1.SZ0
  
  return(VE.S)
}


getVE.Perturb<-function(data,Su,results.Perturb,type,xi){
  nuisance<-results.Perturb$nuisance
  alpha0<-nuisance[[1]]
  alpha1<-nuisance[[2]]
  sigma1<-nuisance[[3]]
  gamma0<-nuisance[[4]]
  gamma1<-nuisance[[5]]
  gamma2<-nuisance[[6]]
  sigma2<-nuisance[[7]]
  muB<-nuisance[[8]]
  muS<-nuisance[[9]]
  sigmaB<-nuisance[[10]]
  sigmaS<-nuisance[[11]]
  rho<-nuisance[[12]]
  sigma3<-sigmaB*sqrt(1-rho^2)
  sigma<-matrix(c(sigmaB^2,rho*sigmaB*sigmaS,rho*sigmaB*sigmaS,sigmaS^2),nrow=2)
  
  data$xi<-xi
  data$S1<-ifelse(data$Z==1,data$S,NA)
  Wu<-sort(unique(data$W[!is.na(data$W)]))
  
  if(type=="marginal"){
    myriskonSW<-riskonSW
    #perturbed version 
    fX<-as.numeric(aggregate(data$xi, by=list(Category=data$W), FUN=sum)[,2])/sum(as.numeric(aggregate(data$xi, by=list(Category=data$W), FUN=sum)[,2]))
    #When s>c
    integrand<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- dnorm(mys,mean=my.mu2,sd=sigma2)*dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    
    fx.s<-NULL
    for(mys in Su[-1]){
      myrow<-NULL
      for(w in 1:4){
        myrow<-c(myrow,integrate(integrand,lower=-Inf,upper=Inf)$value  )
      }
      myrow<-myrow*fX/sum(myrow*fX)
      fx.s<-rbind(fx.s,myrow)
    }
    
    #When s=c
    integrand<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- pnorm(c,mean=my.mu2,sd=sigma2)*dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    myrow<-NULL
    for(w in 1:4){
      myrow<-c(myrow,integrate(integrand,lower=-Inf,upper=Inf)$value  )
    }
    myrow<-myrow*fX/sum(myrow*fX)
    fx.s<-rbind(myrow,fx.s)
  }
  
  if(type=="Bpos"){
    myriskonSW<-riskonSW_Bpos
    #perturbed version 
    #P(W|B>c)
    fX<-as.numeric(aggregate(data$xi, by=list(Category=data$W), FUN=sum)[,2])/sum(as.numeric(aggregate(data$xi, by=list(Category=data$W), FUN=sum)[,2]))
    myrow<-NULL
    for(w in 1:4){
      if(w==1) {my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu1<-alpha0+alpha1[4]}
      myrow<-c(myrow,1-pnorm(c,mean=my.mu1,sd=sigma1))
    }
    fX<-myrow*fX/sum(myrow*fX)
    
    #When s>c
    integrand1<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- dnorm(mys,mean=my.mu2,sd=sigma2)*dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    integrand2<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    
    fx.s<-NULL
    for(mys in Su[-1]){
      myrow<-NULL
      for(w in 1:4){
        myrow<-c(myrow,integrate(integrand1,lower=c,upper=Inf)$value/integrate(integrand2,lower=c,upper=Inf)$value  )
      }
      myrow<-myrow*fX/sum(myrow*fX)
      fx.s<-rbind(fx.s,myrow)
    }
    
    #When s=c
    integrand1<-function(b){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      out<- pnorm(c,mean=my.mu2,sd=sigma2)*dnorm(b,mean=my.mu1,sd=sigma1)
      return(out)
    }
    myrow<-NULL
    for(w in 1:4){
      myrow<-c(myrow,integrate(integrand1,lower=c,upper=Inf)$value/integrate(integrand2,lower=c,upper=Inf)$value  )
    }
    myrow<-myrow*fX/sum(myrow*fX)
    fx.s<-rbind(myrow,fx.s)
  }
  
  if(type=="Bneg"){
    myriskonSW<-riskonSW_Bneg
    #perturbed version 
    #P(W|B=c)
    fX<-as.numeric(aggregate(data$xi, by=list(Category=data$W), FUN=sum)[,2])/sum(as.numeric(aggregate(data$xi, by=list(Category=data$W), FUN=sum)[,2]))
    myrow<-NULL
    for(w in 1:4){
      if(w==1) {my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu1<-alpha0+alpha1[4]}
      myrow<-c(myrow,pnorm(c,mean=my.mu1,sd=sigma1))
    }
    fX<-myrow*fX/sum(myrow*fX)
    
    #When s>c
    b<-c
    fx.s<-NULL
    for(mys in Su[-1]){
      myrow<-NULL
      for(w in 1:4){
        if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
        if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
        if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
        if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
        myrow<-c(myrow, dnorm(mys,mean=my.mu2,sd=sigma2) )
      }
      myrow<-myrow*fX/sum(myrow*fX)
      fx.s<-rbind(fx.s,myrow)
    }
    #When s=c
    b<-c
    mys <-c
    myrow<-NULL
    for(w in 1:4){
      if(w==1) {my.mu2<-gamma0+gamma1*b+gamma2[1]; my.mu1<-alpha0+alpha1[1]}
      if(w==2) {my.mu2<-gamma0+gamma1*b+gamma2[2]; my.mu1<-alpha0+alpha1[2]}
      if(w==3) {my.mu2<-gamma0+gamma1*b+gamma2[3]; my.mu1<-alpha0+alpha1[3]}
      if(w==4) {my.mu2<-gamma0+gamma1*b+gamma2[4]; my.mu1<-alpha0+alpha1[4]}
      myrow<-c(myrow, pnorm(mys,mean=my.mu2,sd=sigma2) )
    }
    myrow<-myrow*fX/sum(myrow*fX)
    fx.s<-rbind(myrow,fx.s)
  }
  
  PY1.SZ0<-PY1.SZ1<-rep(0,length(Su))
  
  for(i in 1:length(Su)){
    for (j in 1:length(Wu)){
      PY1.SZ1[i]=PY1.SZ1[i]+myriskonSW(s=Su[i],w=Wu[j],z=1,nuisance = results.Perturb$nuisance, beta = results.Perturb$beta)*fx.s[i,j]
      PY1.SZ0[i]=PY1.SZ0[i]+myriskonSW(s=Su[i],w=Wu[j],z=0,nuisance = results.Perturb$nuisance, beta = results.Perturb$beta)*fx.s[i,j]
    }
  }
  #in case PY1.SZ0 is very very small, and R recorded as 0:
  PY1.SZ0<-ifelse(PY1.SZ0==0,1e-300,PY1.SZ0)
  VE.S<-1-PY1.SZ1/PY1.SZ0
  return(VE.S)
}


###############################################################################################################################################
#########The following functions are for analysis not incorporating B, to compare efficiency
#Estimate B|W

## (Negative) log-likelihood function
#theta=(alpha0,alpha12,alpha13,alpha14,sigma_1)
logl.SonW <- function(theta,data) {
  #first calculate case&control sampling wright
  w.case <- sum(data$delta_S==1 & data$Y==1 & data$Z==1)/sum(data$Y==1 & data$Z==1)
  w.control <- sum(data$delta_S==1 & data$Y==0 & data$Z==1)/sum(data$Y==0 & data$Z==1)
  
  dat<-data[data$delta_S==1&data$Z==1,]
  alpha0<-theta[1]
  alpha1<-c(0,theta[2:4])
  sigma1<-theta[5]
  mu1<-ifelse(dat$W==1,alpha0+alpha1[1],ifelse(dat$W==2,alpha0+alpha1[2],ifelse(dat$W==3,alpha0+alpha1[3],alpha0+alpha1[4])))
  ff<-ifelse(dat$S==c,pnorm(c,mean=mu1,sd=sigma1), dnorm(dat$S,mean=mu1,sd=sigma1))
  
  wts<- ifelse(dat$Y==1,w.case,w.control)
  wts<-1/wts
  ff<-ff*wts
  return(-sum(log(ff)))
}

###Calculate observed likelihood

like_noB<-function(beta,data,nuisance){
  beta_0<-beta[1];beta_1<-beta[2];beta_2<-beta[3];beta_3<-beta[4];beta_42<-beta[5];beta_43<-beta[6];beta_44<-beta[7]
  
  alpha0<-nuisance[[1]]
  alpha1<-nuisance[[2]]
  sigma1<-nuisance[[3]]
  
  ans <- 0
  
  # delta_S=1,Z=1
  dat<-data[data$delta_S==1&data$Z==1,] #2692
  ff<-ifelse(dat$W==1,beta_0+beta_1+beta_2*dat$S+beta_3*dat$S,
             ifelse(dat$W==2,beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_42,
                    ifelse(dat$W==3, beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_43, beta_0+beta_1+beta_2*dat$S+beta_3*dat$S+beta_44)))
  
  risk<-pnorm(ff)
  rr<-ifelse(dat$Y==1,risk,1-risk)
  ans<-ans+sum(log(rr))
  
  #delta_S=1&Z=0 or delta_S=0
  dat<-data[data$Z==0|data$delta_S==0,] #7308
  mu1<-ifelse(dat$W==1,alpha0+alpha1[1],ifelse(dat$W==2,alpha0+alpha1[2],ifelse(dat$W==3,alpha0+alpha1[3],alpha0+alpha1[4])))
  
  ff<-ifelse(dat$W==1,beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c,
             ifelse(dat$W==2,beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_42,
                    ifelse(dat$W==3, beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_43, beta_0+beta_1*dat$Z+beta_2*c+beta_3*dat$Z*c+beta_44)))
  term1<- pnorm( (c-mu1)/sigma1 ) * pnorm(ff)
  
  myid<-1:nrow(dat)
  integrand<-function(s,id){
    if(dat$W[id]==1) {mu1<-alpha0+alpha1[1]; beta_4<-0}
    if(dat$W[id]==2) {mu1<-alpha0+alpha1[2]; beta_4<-beta_42}
    if(dat$W[id]==3) {mu1<-alpha0+alpha1[3]; beta_4<-beta_43}
    if(dat$W[id]==4) {mu1<-alpha0+alpha1[4]; beta_4<-beta_44}
    
    out<-pnorm(beta_0+beta_1*dat$Z[id]+beta_2*s+beta_3*dat$Z[id]*s+beta_4)*dnorm(s, mean=mu1, sd=sigma1)
    return(out)
  }
  
  area<-function(id) integrate(integrand,lower=c,upper=Inf,id=id)$value
  area2<-Vectorize(area)
  
  term2<-area2(myid)
  
  risk<-term1+term2
  rr<-ifelse(dat$Y==1,risk,1-risk)
  ans<-ans+sum(log(rr))
  
  #(Negative) log-likelihood function
  ans<- -ans
  
}


getBeta_noB<-function(data){
  ################Estimate Nuisance parameters
  #Estimate S|W
  # Estimate theta
  theta.start <- c(3.02,0.23, -0.15, 0.43,0.65)
  out1 <- optim(par=theta.start, fn=logl.SonW,data=data)$par
  #out<- nlm(logl.BonW, theta.start, data = data) produces the same results
  alpha0<-out1[1]
  alpha1<-c(0,out1[2:4])
  sigma1<-out1[5]
  
  nuisance<-list(alpha0,alpha1, sigma1)
  
  # Estimate beta
  beta.start <-c( -0.5442,0.0000, -0.1250, -0.1400,  0.2360,  0.1060,  0.2000)
  out <- optim(par=beta.start, fn=like_noB,data=data,nuisance=nuisance)$par
  out<-list(nuisance=nuisance, beta=out)
  return(out)
}
