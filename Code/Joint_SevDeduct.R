# joint model for claim severity and deductible choice

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
# joint model deductible and severity
# use claim file
#########################################################################

datin.clm$dedratio <- datin.clm$Deduct/datin.clm$Coverage

# response variable
yy1 <- datin.clm$dedratio
yy2 <- datin.clm$LossBeforeDeductible


# read in parameters
gam1 <- readRDS(file="coef.deductGB2.RDS")
gam2 <- readRDS(file="coef.sev.RDS")

# deductible
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage,data=datin.clm))
k1 <- ncol(cova1)

mu1 <- cova1%*%gam1[1:k1]
sigma1 <- gam1[k1+1]
alpha11 <- gam1[k1+2]
alpha21 <- gam1[k1+3]


# severity
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.clm))
k <- ncol(cova) 


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


# now define likelihood for the joint model
# this is two-stage estimation
# severity loglikelihood
loglikSev <- function(parms,fam){
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  rho <- parms[k+4]
  df <- parms[k+5]
  u1 <- pGB2(yy1,mu1,sigma1,alpha11,alpha21)
  u2 <- pGB2(yy2,mu,sigma,alpha1,alpha2)
  llk <- -sum(log(dGB2(yy2,mu,sigma,alpha1,alpha2))) - sum(log(BiCopPDF(u1,u2,fam,rho,df)))  
  llk
}
ini.sev <- c(gam2,0,5)
loglikSev(ini.sev,fam=1)  # this is to check nested case, should be equal to severity marginal

# try different copulas
# fam = 1 gaussian
zop1 <- nlminb(ini.sev,loglikSev,fam=1,lower =c(rep(-Inf,ncol(cova)),-Inf,0.001,0.001,-0.999),upper=c(rep(Inf,ncol(cova)),Inf,Inf,Inf,0.999),control=list(eval.max=500))
print(zop1)
obj <- BiCop(family = 1, par = zop1$par[length(ini.sev)-1], par2 = zop1$par[length(ini.sev)])
summary(obj)


