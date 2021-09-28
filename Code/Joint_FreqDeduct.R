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

datin.pol$dedratio <- datin.pol$Deduct/datin.pol$Coverage

# response variable
yy1 <- datin.pol$dedratio
yy2 <- datin.pol$Freq

# read in parameters
gam1 <- readRDS(file="coef.deductGB2.RDS")
gam2 <- readRDS(file="coef.freq.RDS")


# deductible
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage,data=datin.pol))
k1 <- ncol(cova1)

# define density of GB2
dGB2 <- function(y, mu, sigma, alpha1, alpha2, log = FALSE){
  fstterm <- exp(alpha1 * (log(y) - mu) / sigma)
  sndterm <- y * abs(sigma)* gamma(alpha1) * gamma(alpha2)/gamma(alpha1 + alpha2)
  thdterm <- (1 + exp((log(y) - mu)/sigma ) )^(alpha1 + alpha2)
  ans <- fstterm/(sndterm * thdterm)
  if (log) log(ans) else ans
  ans <- ifelse (is.na(ans), 0, ans)
  ans
}
pGB2 <- function(y, mu, sigma, alpha1, alpha2){
  ndf <- 2*alpha1
  ddf <- 2*alpha2
  r <- (log(y)-mu)/sigma
  z <- (alpha2/alpha1)*exp(r)
  ans <- pf(z,ndf,ddf)
  ans
}

mu <- cova1%*%gam1[1:k1]
sigma <- gam1[k1+1]
alpha1 <- gam1[k1+2]
alpha2 <- gam1[k1+3]


# frequency
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.pol))
k <- ncol(cova)                        # number of predictors

# covariates
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
  u1 <- pGB2(yy1,mu,sigma,alpha1,alpha2)
  u2 <- Fcount(yy2)
  u2a <- Fcount(yy2-1)
  temp <- BiCopHfunc1(u1,u2,fam,rho,df)-BiCopHfunc1(u1,u2a,fam,rho,df)
  llk <- -sum(log(temp))  
  llk
}



# try different copulas
# fam = 1 gaussian
ini.freq <- c(0,5)
zop1 <- nlminb(ini.freq,loglikFreq,fam=1,lower =c(-0.999,1),upper=c(0.999,Inf),control=list(eval.max=500))
print(zop1)
#obj <- BiCop(family = 1, par = zop1$par[length(ini.freq)-1], par2 = zop1$par[length(ini.freq)])
#summary(obj)
#$par -0.2656683  5.0000000
#$objective  3458.754
BiCopPar2Tau(1,-0.2656683)


#readRDS(file="cop.Gaussian.RDS")



#############################################
# The final section is the t copula
#################################################

# fam = 2 t
ini.freq <- c(0,5)
zop <- nlminb(ini.freq,loglikFreq,fam=2,lower =c(-0.999,1),upper=c(0.999,Inf),control=list(eval.max=500))
print(zop)
# $par =  -0.2689659 12.3180148
# loglik =  3453.594 


saveRDS(c(fam=2,zop$par),file="cop.DN.RDS")


