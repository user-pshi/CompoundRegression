# Goodness-of-fit for marginal and conditional models for claim frequency
# Need to run DeductibleModels.R to obtain "coef.deductGB2.RDS"
# Need to run FrequencyModels.R to obtain "coef.freq.RDS"
# need to run Joint_FreqDeduct.R to obtain "cop.DN.RDS"
# Reproduce Table S.1


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

# read in parameters
gam1 <- readRDS(file="coef.deductGB2.RDS")
gam2 <- readRDS(file="coef.freq.RDS")
cop <- readRDS("cop.DN.RDS")



datin.pol$dedratio <- datin.pol$Deduct/datin.pol$Coverage

# response variable
yy1 <- datin.pol$dedratio
yy2 <- datin.pol$Freq



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
k <- ncol(cova)                        

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


# empirical model
numrow <- max(yy2)+1
emp<-rep(0,numrow)
for(i in 1:numrow){
  emp[i]<-sum(yy2==(i-1))
}

# marginal model
exp1nb<-rep(0,numrow)
exp1nb[1]<-sum(pzero+(1-pzero-pone)*dnbinom(0,mu=meannb,size=size))
exp1nb[2]<-sum(pone+(1-pzero-pone)*dnbinom(1,mu=meannb,size=size))
for(i in 3:numrow){
  exp1nb[i]<-sum((1-pzero-pone)*dnbinom(i-1,mu=meannb,size=size))
}

# conditional model
exp1nb.cond<-rep(0,max(yy2)+1)
for(i in 0:max(yy2)){
  u1 <- pGB2(yy1,mu,sigma,alpha1,alpha2)
  u2 <- Fcount(i)
  u2a <- Fcount(i-1)
  exp1nb.cond[i+1] <- sum(BiCopHfunc1(u1,u2,cop[1],cop[2],cop[3])-BiCopHfunc1(u1,u2a,cop[1],cop[2],cop[3]))
}


# reproduce Table S.1
emp2 <- c(emp[1:5],sum(emp[6:7]),sum(emp[8:11]),sum(emp[12:15]),sum(emp[16:numrow]))
exp1nb2 <- c(exp1nb[1:5],sum(exp1nb[6:7]),sum(exp1nb[8:11]),sum(exp1nb[12:15]),sum(exp1nb[16:numrow]))
exp1nb.cond2 <- c(exp1nb.cond[1:5],sum(exp1nb.cond[6:7]),sum(exp1nb.cond[8:11]),sum(exp1nb.cond[12:15]),sum(exp1nb.cond[16:numrow]))
rbind(emp2/sum(emp2)*100,exp1nb2/sum(exp1nb2)*100,exp1nb.cond2/sum(exp1nb.cond2)*100)






