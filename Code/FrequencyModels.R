# Fit models for claim frequency
# Reproduce Table 2


library(ggplot2)
library(mgcv)
library(numDeriv)
library(dglm)
library(ald)
library(quantreg)
library(rcompanion)
library(MASS)

library(BB)
library(MASS)
library(pscl)
library(VineCopula)


rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009; use 2010 as hold-out
selectyear <- 2006:2009
datin.pol <- poldat[which(poldat$Year%in%selectyear),]
datin.clm <- claimdat[which(claimdat$Year%in%selectyear),]



###########################################################################
# Analyze claim frequency 
# use policy file
# Fit 5 candidate regression models
# 1 - Poisson
# 2 - Negative Binomial
# 3 - Zero inflated Poisson
# 4 - Zero inflated Negative Binomial
# 5 - Negative Binomial with inflation in both zero and one
#########################################################################

yy <- datin.pol$Freq
cova <- model.matrix(lm(yy~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.pol))


# covariates
cova.zero <- cova[,c(1,ncol(cova))]
cova.one <- cova[,c(1,ncol(cova))]

n<-ncol(cova)
n0<-ncol(cova.zero)
n1<-ncol(cova.one)



######################################
# Model 1: Possion 
######################################
freq.p <-glm(yy~cova-1,data=datin.pol,family=poisson(link = "log"))
saveRDS(coef(freq.p),file="coef.freq.Poisson.RDS")

######################################
# Model 2: Negative Binomial
######################################
freq.nb<-glm.nb(yy~cova-1,data=datin.pol,link = "log")
saveRDS(c(coef(freq.nb),summary(freq.nb)$theta),file="coef.freq.NB.RDS")

######################################
# Model 3: Zero inflated Poisson
######################################
freq.zerop<-zeroinfl(yy~ cova-1|cova.zero-1, data=datin.pol,dist = c("poisson"),link = c("logit"))


######################################
# Model 4: Zero inflated Negative Binomial
######################################
freq.zeronb<-zeroinfl(yy~ cova-1|cova.zero-1, data=datin.pol,dist = c("negbin"),link = c("logit"))


##################################################
# Model 5: Negative Binomial with inflation in both zero and one
#################################################
loglik01nb<-function(a){
  beta<-a[1:n]
  gamma0<-a[(n+1):(n+n0)]
  gamma1<-a[(n+n0+1):(n+n0+n1)]
  size<-a[n+n0+n1+1]
  reg0<-cova.zero%*%gamma0
  reg1<-cova.one%*%gamma1
  pzero<-1/(exp(-reg0)+exp(reg1-reg0)+1)
  pone<-1/(exp(-reg1)+exp(reg0-reg1)+1)
  reg<-pmin(cova%*%beta,23)
  meannb<-exp(reg)
  l<--sum(log(pmax(pzero*(yy==0)+pone*(yy==1)+(1-pzero-pone)*dnbinom(yy,mu=meannb,size=size))))  
  return(l)
}

initial <- c(freq.zeronb$coefficients$count,freq.zeronb$coefficients$zero,rep(0,n1),summary(freq.zeronb)$theta)
zop <- nlminb(initial,loglik01nb)  
zop
hess<-hessian(loglik01nb,zop$par)
se <-sqrt(diag(solve(hess)))

################################################################
# Estimate and save parameter estimates for frequency component
# reproduce Table 2: frequency component
##############################################################
cbind(zop$par,se)


saveRDS(zop$par,file="coef.freq.RDS")




