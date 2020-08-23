#set true parameter values 
alpha0<-1.728211
alpha1<-c(-0.35, 0.58, 0.9, -0.6)
sigma1<-0.86

gamma0<-1.5
gamma1<-0.5
gamma2<-c(0,0.2,-0.1,0.4)
sigma2<-0.4

muB<-alpha0+alpha1
muS<-c(gamma0+gamma1*(alpha0+alpha1[1])+gamma2[1],
       gamma0+gamma1*(alpha0+alpha1[2])+gamma2[2],
       gamma0+gamma1*(alpha0+alpha1[3])+gamma2[3],
       gamma0+gamma1*(alpha0+alpha1[4])+gamma2[4])
sigmaB<-sigma1
sigmaS<- sqrt(sigma2^2+ (gamma1^2)*(sigma1^2))
rho<- gamma1*sigma1/sqrt(sigma2^2+gamma1^2*sigma1^2)

nuisance_true<-list(alpha0,alpha1, sigma1, gamma0, gamma1, gamma2, sigma2, muB, muS, sigmaB, sigmaS, rho)

beta_true<-c(-0.10, 0.16,  -0.34 , -0.21, -0.25, 0.2360,0.1060,0.2000)
c<-1