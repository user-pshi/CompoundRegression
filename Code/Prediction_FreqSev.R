# Hold-out sample prediction for claim frequency and severity
# Need to run DeductibleModels.R to obtain "coef.deductGB2.RDS"
# Need to run FrequencyModels.R to obtain "coef.freq.RDS"
# Need to run Joint_FreqDeduct.R to obtain "cop.DN.RDS"
# Need to run JointModel.R to obtain "cop.DS.RDS"
# reproduce results in Table 3


library(ggplot2)
library(mgcv)
library(numDeriv)
library(dglm)
library(ald)
library(quantreg)
library(rcompanion)
library(MASS)
library(cplm)


library(BB)
library(MASS)
library(pscl)
library(VineCopula)


rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009; use 2010 as hold-out
selectyear <- 2010
datout.pol <- poldat[which(poldat$Year%in%selectyear),]
datout.clm <- claimdat[which(claimdat$Year%in%selectyear),]


# read in parameters
gam1 <- readRDS(file="coef.deductGB2.RDS")
gam2 <- readRDS(file="coef.freq.RDS")
copDN <- readRDS(file="cop.DN.RDS")
copDS <- readRDS(file="cop.DS.RDS")




##################################################
# calculate score for frequency model
##################################################

datout.pol$dedratio <- datout.pol$Deduct/datout.pol$Coverage

# response variable
yy1 <- datout.pol$dedratio
yy2 <- datout.pol$Freq



# deductible
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage,data=datout.pol))
k1 <- ncol(cova1)


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

cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datout.pol))
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



# cdf of 01 inflated NB
Fcount <- function(y){
  ans <- pzero*(y>=0)+pone*(y>=1)+(1-pzero-pone)*pnbinom(y,mu=meannb,size=size)
  ans
}



# independence density
fcount.ind <- function(y){
  ans <- pzero*(y==0)+pone*(y==1)+(1-pzero-pone)*dnbinom(y,mu=meannb,size=size)
  ans
}

# dependence density
fcount.dep <- function(y){
  u1 <- pGB2(yy1,mu,sigma,alpha1,alpha2)
  u2 <- Fcount(y)
  u2a <- Fcount(y-1)
  ans <- BiCopHfunc1(u1,u2,family=copDN[1],par=copDN[2],par2=copDN[3])-BiCopHfunc1(u1,u2a,family=copDN[1],par=copDN[2],par2=copDN[3])
  ans
}

fscore.ind <- fcount.ind(yy2)
fscore.dep <- fcount.dep(yy2)
f.ind <- sum(log(fscore.ind))
f.dep <- sum(log(fscore.dep))
f.perc <- sum((fscore.dep>=fscore.ind))/length(fscore.dep)



#################################################
# calculate score for severity model
##################################################

datout.clm$dedratio <- datout.clm$Deduct/datout.clm$Coverage

# response variable
yy1 <- datout.clm$dedratio
yy2 <- datout.clm$Freq
yy3 <- datout.clm$LossBeforeDeductible

# deductible
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage,data=datout.clm))
k1 <- ncol(cova1)

mu1 <- cova1%*%gam1[1:k1]
sigma1 <- gam1[k1+1]
alpha11 <- gam1[k1+2]
alpha21 <- gam1[k1+3]


# covariatesfor frequency and severity
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datout.clm))
k <- ncol(cova) 

# frequency
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
  ans
}
# cdf of 01 inflated NB
Fcount <- function(y){
  ans <- pzero*(y>=0)+pone*(y>=1)+(1-pzero-pone)*pnbinom(y,mu=meannb,size=size)
  ans
}


# severity

# define density of GB2
dGB2 <- function(y, mu, sigma, alpha1, alpha2, log = FALSE){
  fstterm <- exp(alpha1 * (log(y) - mu) / sigma)
  sndterm <- y * abs(sigma)* gamma(alpha1) * gamma(alpha2)/gamma(alpha1 + alpha2)
  thdterm <- (1 + exp((log(y) - mu)/sigma ))^(alpha1 + alpha2)
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


famount.ind <- function(y){
  mu <- cova%*%copDS[3:(k+2)]
  sigma <- copDS[k+3]
  alpha1 <- copDS[k+4]
  alpha2 <- copDS[k+5]
  ans <- dGB2(y,mu,sigma,alpha1,alpha2)
  ans 
}


famount.dep <- function(y){
  fam <- copDS[1:2]
  mu <- cova%*%copDS[3:(k+2)]
  sigma <- copDS[k+3]
  alpha1 <- copDS[k+4]
  alpha2 <- copDS[k+5]
  rho1 <- copDS[k+6]
  rho2 <- copDS[k+7]
  df1 <- 10
  df2 <- 10
  
  # frequency|deductible
  u1 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2),family=copDN[1],par=copDN[2],par2=copDN[3])
  u1a <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2-1),family=copDN[1],par=copDN[2],par2=copDN[3])
  # severity|deductible
  u2 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),pGB2(y,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
  # density severity|deductible
  temp1 <- dGB2(y,mu,sigma,alpha1,alpha2)
  temp2 <- BiCopPDF(pGB2(yy1,mu1,sigma1,alpha11,alpha21),pGB2(y,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
  # (frequency,severity)|deductible
  temp3 <- BiCopHfunc2(u1,u2,fam[2],rho2,df2)-BiCopHfunc2(u1a,u2,fam[2],rho2,df2)
  
  ans <- temp1*temp2*temp3/(u1-u1a)
  ans
}

sscore.ind <- famount.ind(yy3)
sscore.dep <- famount.dep(yy3)
s.ind <- sum(log(sscore.ind))
s.dep <- sum(log(sscore.dep))
s.perc <- sum((sscore.dep>=sscore.ind))/length(sscore.dep)


# reproduce Table 3
cbind(percent = c(f.perc,s.perc), exogenous = c(f.ind,s.ind), endogenous = c(f.dep,s.dep))
