# Fit models for claim frequency

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
#########################################################################

yy <- datin.pol$Freq
cova <- model.matrix(lm(yy~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.pol))


# covariates
cova.zero <- cova[,c(1,ncol(cova))]
cova.one <- cova[,c(1,ncol(cova))]

n<-ncol(cova)
n0<-ncol(cova.zero)
n1<-ncol(cova.one)


# empirical model
numrow<-max(yy)+1
emp<-rep(0,numrow)
for(i in 1:numrow){
  emp[i]<-sum(yy==(i-1))
}

# possion for comparision
freq.p <-glm(yy~cova-1,data=datin.pol,family=poisson(link = "log"))
beta<-freq.p$coefficients
meanpoisson<-freq.p$fitted.values
expp<-rep(0,numrow)
for(i in 1:numrow){
  expp[i]<-sum(dpois((i-1),meanpoisson))
}


saveRDS(coef(freq.p),file="coef.freq.Poisson.RDS")




# nb for comparision
freq.nb<-glm.nb(yy~cova-1,data=datin.pol,link = "log")
beta<-freq.nb$coefficients
meannb<-freq.nb$fitted.values
expnb<-rep(0,numrow)
for(i in 1:numrow){
  expnb[i]<-sum(dnbinom(i-1,mu=meannb,size=summary(freq.nb)$theta))
}


saveRDS(c(coef(freq.nb),summary(freq.nb)$theta),file="coef.freq.NB.RDS")


# zeroinf poisson
freq.zerop<-zeroinfl(yy~ cova-1|cova.zero-1, data=datin.pol,dist = c("poisson"),link = c("logit"))
gamma<-freq.zerop$coefficients$zero
beta<-freq.zerop$coefficients$count
pzero<-exp(cova.zero%*%gamma)/(1+exp(cova.zero%*%gamma))
meanpoisson<-exp(cova%*%beta)
exp0p<-rep(0,numrow)
exp0p[1]<-sum(pzero+(1-pzero)*dpois(0,meanpoisson))
for(i in 2:numrow){
  exp0p[i]<-sum((1-pzero)*dpois((i-1),meanpoisson))
}

# zeroinflated negative binomial 
freq.zeronb<-zeroinfl(yy~ cova-1|cova.zero-1, data=datin.pol,dist = c("negbin"),link = c("logit"))
gamma<-freq.zeronb$coefficients$zero
beta<-freq.zeronb$coefficients$count
pzero<-exp(cova.zero%*%gamma)/(1+exp(cova.zero%*%gamma))
meannb<-exp(cova%*%beta)
exp0nb<-rep(0,numrow)
exp0nb[1]<-sum(pzero+(1-pzero)*dnbinom(0,mu=meannb,size=freq.zeronb$theta))
for(i in 2:numrow){
  exp0nb[i]<-sum((1-pzero)*dnbinom(i-1,mu=meannb,size=freq.zeronb$theta))
}



# 01 inflated poisson
loglik01p<-function(a){
  beta<-a[1:n]
  gamma0<-a[(n+1):(n+n0)]
  gamma1<-a[(n+n0+1):(n+n0+n1)]
  reg0<-cova.zero%*%gamma0
  reg1<-cova.one%*%gamma1
  pzero<-1/(exp(-reg0)+exp(reg1-reg0)+1)
  pone<-1/(exp(-reg1)+exp(reg0-reg1)+1)
  reg<-cova%*%beta
  meanpoisson<-exp(reg)
  l<--sum(log(pmax((pzero)*(yy==0)+pone*(yy==1)+(1-pzero-pone)*dpois(yy,meanpoisson),10^-10)))
  return(l)
}

initial <- c(freq.zerop$coefficients$count,freq.zerop$coefficients$zero,rep(0,n1))
zop <- nlminb(initial,loglik01p)  
zop

beta<-zop$par[1:n]
gamma0<-zop$par[(n+1):(n+n0)]
gamma1<-zop$par[(n+n0+1):(n+n0+n1)]
pzero<-exp(cova.zero%*%gamma0)/(1+exp(cova.zero%*%gamma0)+exp(cova.one%*%gamma1))
pone<-exp(cova.one%*%gamma1)/(1+exp(cova.one%*%gamma1)+exp(cova.zero%*%gamma0))
meanpoisson<-exp(cova%*%beta)
exp1p<-rep(0,numrow)
exp1p[1]<-sum(pzero+(1-pzero-pone)*dpois(0,meanpoisson))
exp1p[2]<-sum(pone+(1-pzero-pone)*dpois(1,meanpoisson))
for(i in 3:numrow){
  exp1p[i]<-sum((1-pzero-pone)*dpois((i-1),meanpoisson))
}


# 01 inflated nb
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
  l<--sum(log(pmax(pzero*(yy==0)+pone*(yy==1)+(1-pzero-pone)*dnbinom(yy,mu=meannb,size=size),10^-4)))  
  return(l)
}

initial <- c(freq.zeronb$coefficients$count,freq.zeronb$coefficients$zero,rep(0,n1),summary(freq.zeronb)$theta)
zop <- nlminb(initial,loglik01nb)  
zop


saveRDS(zop$par,file="coef.freq.RDS")


#zop$objective
#loglik =  3542.849

# get standard error
library(numDeriv)
hess<-hessian(loglik01nb,zop$par)
se <-sqrt(diag(solve(hess)))
tratio <- zop$par/se
cbind(zop$par,se,tratio)


beta<-zop$par[1:n]
gamma0<-zop$par[(n+1):(n+n0)]
gamma1<-zop$par[(n+n0+1):(n+n0+n1)]
#logsize=zop$par[n+n0+n1+1]
size=zop$par[n+n0+n1+1]
pzero<-exp(cova.zero%*%gamma0)/(1+exp(cova.zero%*%gamma0)+exp(cova.one%*%gamma1))
pone<-exp(cova.one%*%gamma1)/(1+exp(cova.one%*%gamma1)+exp(cova.zero%*%gamma0))
meannb<-exp(cova%*%beta)
exp1nb<-rep(0,numrow)
exp1nb[1]<-sum(pzero+(1-pzero-pone)*dnbinom(0,mu=meannb,size=size))
exp1nb[2]<-sum(pone+(1-pzero-pone)*dnbinom(1,mu=meannb,size=size))
for(i in 3:numrow){
  exp1nb[i]<-sum((1-pzero-pone)*dnbinom(i-1,mu=meannb,size=size))
}






#expected frequency table
table(datin.pol$Freq)
cutoff <- 9
round(cbind(c(emp[1:cutoff],sum(emp[(cutoff+1):numrow])),            
            c(expp[1:cutoff],sum(expp[(cutoff+1):numrow])),
            c(expnb[1:cutoff],sum(expnb[(cutoff+1):numrow])),
            c(exp0p[1:cutoff],sum(exp0p[(cutoff+1):numrow])),
            c(exp0nb[1:cutoff],sum(exp0nb[(cutoff+1):numrow])),
            c(exp1p[1:cutoff],sum(exp1p[(cutoff+1):numrow])),
            c(exp1nb[1:cutoff],sum(exp1nb[(cutoff+1):numrow]))),digits=3)

# percentage of zero
round(cbind(emp,expp,expnb,exp0p,exp0nb,exp1p,exp1nb)[1,]/
        colSums(cbind(emp,expp,expnb,exp0p,exp0nb,exp1p,exp1nb)),digits=3)

# percentage of one
round(cbind(emp,expp,expnb,exp0p,exp0nb,exp1p,exp1nb)[2,]/
        colSums(cbind(emp,expp,expnb,exp0p,exp0nb,exp1p,exp1nb)),digits=3)

#goodness of fit test
round(colSums((cbind(c(expp[1:cutoff],sum(expp[(cutoff+1):numrow])),
                     c(expnb[1:cutoff],sum(expnb[(cutoff+1):numrow])),
                     c(exp0p[1:cutoff],sum(exp0p[(cutoff+1):numrow])),
                     c(exp0nb[1:cutoff],sum(exp0nb[(cutoff+1):numrow])),
                     c(exp1p[1:cutoff],sum(exp1p[(cutoff+1):numrow])),
                     c(exp1nb[1:cutoff],sum(exp1nb[(cutoff+1):numrow])))
               -c(emp[1:cutoff],sum(emp[(cutoff+1):numrow])))^2/
                (cbind(
                  c(expp[1:cutoff],sum(expp[(cutoff+1):numrow])),
                  c(expnb[1:cutoff],sum(expnb[(cutoff+1):numrow])),
                  c(exp0p[1:cutoff],sum(exp0p[(cutoff+1):numrow])),
                  c(exp0nb[1:cutoff],sum(exp0nb[(cutoff+1):numrow])),
                  c(exp1p[1:cutoff],sum(exp1p[(cutoff+1):numrow])),
                  c(exp1nb[1:cutoff],sum(exp1nb[(cutoff+1):numrow]))))),digits=3)


# likelihood function



