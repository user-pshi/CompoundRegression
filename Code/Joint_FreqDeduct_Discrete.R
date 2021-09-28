# joint model for claim frequency and deductible choice


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
rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009; use 2010 as hold-out
selectyear <- 2006:2009
datin.pol <- poldat[which(poldat$Year%in%selectyear),]
datin.clm <- claimdat[which(claimdat$Year%in%selectyear),]


###########################################################################
# joint model deductible and frequency
# use policy file
#########################################################################

names(datin.pol)
table(datin.pol$Deduct)
datin.pol$deduct <- 1*(datin.pol$Deduct==500) + 2*(datin.pol$Deduct==1000) + 3*(datin.pol$Deduct==2500) +
  4*(datin.pol$Deduct==5000) + 5*(datin.pol$Deduct==10000) + 6*(datin.pol$Deduct==15000) +
  7*(datin.pol$Deduct==25000) + 8*(datin.pol$Deduct==50000) + 9*(datin.pol$Deduct==75000) +
  10*(datin.pol$Deduct==100000)
table(datin.pol$deduct)

# response variable
yy1 <- datin.pol$deduct
yy2 <- datin.pol$Freq

# read in parameters
gam1 <- readRDS(file="coef.deductClogit.RDS")
gam2 <- readRDS(file="coef.freq.RDS")


                      # number of predictors

######################
# deductible
######################
# covariates
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage)-1,data=datin.pol))

a <- gam1[1:9]
b <- gam1[10:15]
p1 <- exp(a[1] + cova1%*%b)/(1+exp(a[1] + cova1%*%b))
p2 <- exp(a[2] + cova1%*%b)/(1+exp(a[2] + cova1%*%b))
p3 <- exp(a[3] + cova1%*%b)/(1+exp(a[3] + cova1%*%b))
p4 <- exp(a[4] + cova1%*%b)/(1+exp(a[4] + cova1%*%b))
p5 <- exp(a[5] + cova1%*%b)/(1+exp(a[5] + cova1%*%b))
p6 <- exp(a[6] + cova1%*%b)/(1+exp(a[6] + cova1%*%b))
p7 <- exp(a[7] + cova1%*%b)/(1+exp(a[7] + cova1%*%b))
p8 <- exp(a[8] + cova1%*%b)/(1+exp(a[8] + cova1%*%b))
p9 <- exp(a[9] + cova1%*%b)/(1+exp(a[9] + cova1%*%b))

# pdf of cumulative logit
fclogit <- function(y){
  mp <- (y==1)*p1 + (y==2)*(p2-p1) + (y==3)*(p3-p2) + (y==4)*(p4-p3) + (y==5)*(p5-p4) + (y==6)*(p6-p5) +
    (y==7)*(p7-p6)  + (y==8)*(p8-p7) +  (y==9)*(p9-p8) + (y==10)*(1-p9)
  ans <- ifelse(y<1,0,mp)
  ans
}
# cdf of cumulative logit
Fclogit <- function(y){
  cp <- (y==1)*p1 + (y==2)*p2 + (y==3)*p3 + (y==4)*p4 + (y==5)*p5 + (y==6)*p6 +
         (y==7)*p7 + (y==8)*p8 + (y==9)*p9 + (y==10)*1
  ans <- ifelse(y<1,0,cp)
  ans
}


##################
# frequency
####################
# covariates
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.pol))

cova.zero <- cova[,c(1,ncol(cova))]
cova.one <- cova[,c(1,ncol(cova))]

n<-ncol(cova)
n0<-ncol(cova.zero)
n1<-ncol(cova.one)

beta<-gam2[1:n]
gamma0<-gam2[(n+1):(n+n0)]
gamma1<-gam2[(n+n0+1):(n+n0+n1)]
size<-gam2[n+n0+n1+1]

pzero<-exp(cova.zero%*%gamma0)/(1+exp(cova.zero%*%gamma0)+exp(cova.one%*%gamma1))
pone<-exp(cova.one%*%gamma1)/(1+exp(cova.one%*%gamma1)+exp(cova.zero%*%gamma0))
meannb<-exp(cova%*%beta)

# pdf of 01 inflated NB
fcount <- function(y){
  ans <- pzero*(y==0)+pone*(y==1)+(1-pzero-pone)*dnbinom(y,mu=meannb,size=size)
  #ans <- pmax(1e-10,ans)
  ans
}
# cdf of 01 inflated NB
Fcount <- function(y){
  ans <- pzero*(y>=0)+pone*(y>=1)+(1-pzero-pone)*pnbinom(y,mu=meannb,size=size)
  ans
}


# now define likelihood for the joint model

loglikFreq <- function(parms,fam){
  rho <- parms[1]
  df <- parms[2]
  u1 <- Fclogit(yy1)
  u1a <- Fclogit(yy1-1)
  u2 <- Fcount(yy2)
  u2a <- Fcount(yy2-1)
  temp <- BiCopCDF(u1,u2,fam,rho,df)-BiCopCDF(u1a,u2,fam,rho,df)-BiCopCDF(u1,u2a,fam,rho,df)+BiCopCDF(u1a,u2a,fam,rho,df)
  llk <- -sum(log(temp))  
  llk
}





#############################################
# The final section is the t copula
#################################################

# fam = 24 Gumbel 90 
ini.freq <- c(-1.5,5)
zop <- nlminb(ini.freq,loglikFreq,fam=24,lower =c(-100,1),upper=c(-1.001,Inf),control=list(eval.max=500))
print(zop)
#$par  -1.278045
#$objective  8806.387
BiCopPar2Tau(24, -1.278045)
#tau  -0.2175549


# get standard error
library(numDeriv)
hess<-hessian(loglikFreq,zop$par,fam=24)
se <-sqrt(diag(solve(hess[1,1])))
tratio <- zop$par[1]/se
cbind(zop$par[1],se,tratio)
# se    tratio
#[1,] -1.278045 0.0224594 -56.90471


saveRDS(c(fam=24,zop$par),file="cop.DN_Discrete.RDS")


