# Gaussian copula model in robust analysis
# Freq: Negative Binomial + Severity: LogNormal
# Need to run DeductibleModels.R to obtain "coef.deductGB2.RDS"
# Need to run FrequencyModels.R to obtain "coef.freq.NB.RDS"
# reproduce copula estimates in Table S.5

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
gam2 <- readRDS(file="coef.freq.NB.RDS")


###############################################################
# Estimate the joint model for deductible and claim frequency
# use policy file
#############################################################

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
n<-ncol(cova)
beta<-gam2[1:n]
size<-gam2[n+1]

meannb<-exp(cova%*%beta)

# pdf of NegBin
fcount <- function(y){
  ans <- dnbinom(y,mu=meannb,size=size)
  #ans <- pmax(1e-10,ans)
  ans
}
# cdf of NegBin
Fcount <- function(y){
  ans <- pnbinom(y,mu=meannb,size=size)
  ans
}


# define likelihood for the joint model

loglikFreq <- function(parms,fam){
  rho <- parms[1]
  df <- parms[2]
  u1 <- pGB2(yy1,mu,sigma,alpha1,alpha2)
  u2 <- Fcount(yy2)
  u2a <- Fcount(yy2-1)
  temp <- BiCopHfunc1(u1,u2,fam,rho,df)-BiCopHfunc1(u1,u2a,fam,rho,df)
  temp <- temp 
  llk <- -sum(log(temp))  
  llk
}



ini.freq <- c(0,5)
zop1 <- nlminb(ini.freq,loglikFreq,fam=1,lower =c(-0.999,1),upper=c(0.999,Inf),control=list(eval.max=500))
print(zop1)

hess1<-hessian(loglikFreq,zop1$par,fam=1)
se1 <-sqrt(diag(solve(hess1[1,1])))


# define Gaussian copula
copDN <- c(1,zop1$par)




###########################################################################
# joint model of frequency and severity conditional on deductible
# use claim file
#########################################################################

datin.clm$dedratio <- datin.clm$Deduct/datin.clm$Coverage

# response variable
yy1 <- datin.clm$dedratio
yy2 <- datin.clm$Freq
yy3 <- datin.clm$LossBeforeDeductible


# deductible
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage,data=datin.clm))
k1 <- ncol(cova1)

mu1 <- cova1%*%gam1[1:k1]
sigma1 <- gam1[k1+1]
alpha11 <- gam1[k1+2]
alpha21 <- gam1[k1+3]


# covariates for frequency and severity
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.clm))
k <- ncol(cova) 

# frequency
# covariates
n<-ncol(cova)
beta<-gam2[1:n]
size<-gam2[n+1]

meannb<-exp(cova%*%beta)

# pdf of NegBin
fcount <- function(y){
  ans <- dnbinom(y,mu=meannb,size=size)
  #ans <- pmax(1e-10,ans)
  ans
}
# cdf of NegBin
Fcount <- function(y){
  ans <- pnbinom(y,mu=meannb,size=size)
  ans
}



# severity
# get an initial estimate
# LN
fit1 <- lm(log(yy3)~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.clm)
gam3 <- c(coef(fit1),summary(fit1)$sigma)



# define likelihood for the joint model
loglikJoint1 <- function(parms,fam){
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  rho1 <- parms[k+2]
  rho2 <- parms[k+3]
  df1 <- 5 
  df2 <- 5 
  
  # frequency|deductible
  u1 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2),family=copDN[1],par=copDN[2],par2=copDN[3])
  u1a <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2-1),family=copDN[1],par=copDN[2],par2=copDN[3])
  # severity|deductible
  u2 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),plnorm(yy3,meanlog=mu,sdlog=sigma),fam[1],rho1,df1)
  # density severity|deductible
  temp1 <- dlnorm(yy3,meanlog=mu,sdlog=sigma)
  temp2 <- BiCopPDF(pGB2(yy1,mu1,sigma1,alpha11,alpha21),plnorm(yy3,meanlog=mu,sdlog=sigma),fam[1],rho1,df1)
  # (frequency,severity)|deductible
  temp3 <- BiCopHfunc2(u1,u2,fam[2],rho2,df2)-BiCopHfunc2(u1a,u2,fam[2],rho2,df2)
  
  # - log likelihood
  llk <- -sum(log(temp1)) - sum(log(temp2)) - sum(log(temp3))
  llk
  
  # residual of conditional severity
  resid <- qnorm((BiCopCDF(u1,u2,fam[2],rho2,df2) - BiCopCDF(u1a,u2,fam[2],rho2,df2))/(u1-u1a))
  
  ans <- list(llk=llk,resid=resid)
  ans
  
}
loglikJoint <- function(parms,fam) loglikJoint1(parms,fam)$llk



ini.joint <- c(gam3,0,0)
zop <- nlminb(ini.joint,loglikJoint,fam=c(1,1),lower =c(rep(-Inf,ncol(cova)),0.001,-0.999,-0.999),upper=c(rep(Inf,ncol(cova)),Inf,0.999,0.999),control=list(eval.max=500))
print(zop)
hess<-hessian(loglikJoint,zop$par,fam=c(1,1))
se <-sqrt(diag(solve(hess)))



# estimated copula parameter in Table S.5
cbind(c(zop1$par[1],zop$par[9:10]),c(se1,se[9:10]))






