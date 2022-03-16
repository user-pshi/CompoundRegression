# Compare simplified copula and conditional copula models
# Need to run DeductibleModels.R to obtain "coef.deductGB2.RDS"
# Need to run FrequencyModels.R to obtain "coef.freq.RDS"
# reproduce results in Table S.4

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


###########################################################################
# joint model: deductible and frequency
# 
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


#############################################
# Effect of Entity Type
#############################################

loglikFreq <- function(parms,fam){
  df <- parms[1]
  z <- parms[2] + parms[3]*((datin.pol$TypeCity==1) + (datin.pol$TypeCounty==1) + (datin.pol$TypeTown==1) + (datin.pol$TypeVillage==1))
  rho <- (exp(2*z)-1)/(exp(2*z)+1)
  u1 <- pGB2(yy1,mu,sigma,alpha1,alpha2)
  u2 <- Fcount(yy2)
  u2a <- Fcount(yy2-1)
  temp <- BiCopHfunc1(u1,u2,fam,rho,df)-BiCopHfunc1(u1,u2a,fam,rho,df)
  llk <- -sum(log(temp))  
  llk
}

# fam = 2 t
ini.freq <- c(5,1,1)
zop <- nlminb(ini.freq,loglikFreq,fam=2, lower =c(1,-Inf,-Inf),upper=c(Inf,Inf,Inf),control=list(eval.max=500))
print(zop)

l1.entity <- zop$objective
parms.entity <- zop$par



#############################################
# Effect of Coverage
#############################################
loglikFreq <- function(parms,fam){
  df <- parms[1]
  z <- parms[2] + parms[3]*datin.pol$Coverage/1000
  rho <- (exp(2*z)-1)/(exp(2*z)+1)
  u1 <- pGB2(yy1,mu,sigma,alpha1,alpha2)
  u2 <- Fcount(yy2)
  u2a <- Fcount(yy2-1)
  temp <- BiCopHfunc1(u1,u2,fam,rho,df)-BiCopHfunc1(u1,u2a,fam,rho,df)
  llk <- -sum(log(temp))  
  llk
}


# fam = 2 t
ini.freq <- c(5,1,1,0)
zop <- nlminb(ini.freq,loglikFreq,fam=2, lower =c(1,-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf),control=list(eval.max=500))
print(zop)

l1.coverage <- zop$objective
parms.coverage <- zop$par

#############################################
# Effect of Entity type + Coverage
#############################################

loglikFreq <- function(parms,fam){
  df <- parms[1]
  z <- parms[2] + parms[3]*datin.pol$Coverage/1000 + 
    parms[4]*((datin.pol$TypeCity==1) + (datin.pol$TypeCounty==1) + (datin.pol$TypeTown==1) + (datin.pol$TypeVillage==1))
  rho <- (exp(2*z)-1)/(exp(2*z)+1)
  u1 <- pGB2(yy1,mu,sigma,alpha1,alpha2)
  u2 <- Fcount(yy2)
  u2a <- Fcount(yy2-1)
  temp <- BiCopHfunc1(u1,u2,fam,rho,df)-BiCopHfunc1(u1,u2a,fam,rho,df)
  llk <- -sum(log(temp))  
  llk
}

# fam = 2 t
ini.freq <- c(5,1,1,0,0)
zop <- nlminb(ini.freq,loglikFreq,fam=2, lower =c(1,-Inf,-Inf,-Inf,-Inf),upper=c(Inf,Inf,Inf,Inf,Inf),control=list(eval.max=500))
print(zop)

l1.both <- zop$objective
parms.both <- zop$par




###########################################################################
# joint model of frequency and severity conditional on deductible
# 
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


# covariatesfor frequency and severity
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.clm))
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

# get an initial estimate
# define loglikelihood function
loglikM <- function(parms){         
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  llk <- -sum(log(dGB2(yy3, mu, sigma, alpha1, alpha2)))
  llk
}
fit <- lm(log(yy3)~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.clm)
init <- c(fit$coefficients,1,1,1)
zop <- nlminb(init,loglikM, lower =c(rep(-Inf,ncol(cova)),-Inf,0.001,0.001))   
gam3 <- zop$par



#############################################
# Effect of Entity Type
#############################################

# Define association parameter in copula (R,N)
z <- parms.entity[2] + parms.entity[3]*((datin.clm$TypeCity==1) + (datin.clm$TypeCounty==1) + (datin.clm$TypeTown==1) + (datin.clm$TypeVillage==1))
rho.entity <- (exp(2*z)-1)/(exp(2*z)+1)

# Define likelihood for the joint model
loglikJoint1 <- function(parms,fam){
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  df1 <- 5
  df2 <- 5
  
  z1 <- parms[k+4] + parms[k+5]*((datin.clm$TypeCity==1) + (datin.clm$TypeCounty==1) + (datin.clm$TypeTown==1) + (datin.clm$TypeVillage==1))
  z2 <- parms[k+6] + parms[k+7]*((datin.clm$TypeCity==1) + (datin.clm$TypeCounty==1) + (datin.clm$TypeTown==1) + (datin.clm$TypeVillage==1))
  rho1 <- z1
  rho2 <- -exp(z2)
  
  # frequency|deductible
  u1 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2),family=2,par=rho.entity,par2=parms.entity[1])
  u1a <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2-1),family=2,par=rho.entity,par2=parms.entity[1])
  # severity|deductible
  u2 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
  # density severity|deductible
  temp1 <- dGB2(yy3,mu,sigma,alpha1,alpha2)
  temp2 <- BiCopPDF(pGB2(yy1,mu1,sigma1,alpha11,alpha21),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
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



# The final model is cop1 = Frank, cop2 = rotated Gubmel 
ini.joint <- c(gam3,1,0,1,0)
loglikJoint(ini.joint,fam=c(5,24))  #31703.93
zop <- nlminb(ini.joint,loglikJoint,fam=c(5,24),control=list(eval.max=500))
print(zop)

l2.entity <- zop$objective
l.entity <- l1.entity + l2.entity



#############################################
# Effect of Coverage
#############################################

# Define association parameter in copula (R,N)
z <- parms.coverage[2] + parms.coverage[3]*datin.clm$Coverage/1000
rho.coverage <- (exp(2*z)-1)/(exp(2*z)+1)

# Define likelihood for the joint model
loglikJoint1 <- function(parms,fam){
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  df1 <- 5
  df2 <- 5
  
  z1 <- parms[k+4] + parms[k+5]*datin.clm$Coverage/1000
  z2 <- parms[k+6] + parms[k+7]*datin.clm$Coverage/1000
  rho1 <- z1
  rho2 <- -exp(z2)
  
  # frequency|deductible
  u1 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2),family=2,par=rho.coverage,par2=parms.coverage[1])
  u1a <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2-1),family=2,par=rho.coverage,par2=parms.coverage[1])
  # severity|deductible
  u2 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
  # density severity|deductible
  temp1 <- dGB2(yy3,mu,sigma,alpha1,alpha2)
  temp2 <- BiCopPDF(pGB2(yy1,mu1,sigma1,alpha11,alpha21),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
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



# The final model is cop1 = Frank, cop2 = rotated Gubmel 
ini.joint <- c(gam3,1,0,1,0)
loglikJoint(ini.joint,fam=c(5,24))  #31708.51
zop <- nlminb(ini.joint,loglikJoint,fam=c(5,24),control=list(eval.max=500))
print(zop)

l2.coverage <- zop$objective
l.coverage <- l1.coverage + l2.coverage


#############################################
# Effect of Entity type + Coverage
#############################################

# Define association parameter in copula (R,N)
z <- parms.both[2] + parms.both[3]*datin.clm$Coverage/1000 + 
     parms.both[4]*((datin.clm$TypeCity==1) + (datin.clm$TypeCounty==1) + (datin.clm$TypeTown==1) + (datin.clm$TypeVillage==1))
rho.both <- (exp(2*z)-1)/(exp(2*z)+1)

# Define likelihood for the joint model
loglikJoint1 <- function(parms,fam){
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  df1 <- 5
  df2 <- 5
  
    z1 <- parms[k+4] + parms[k+5]*datin.clm$Coverage/1000 + 
          parms[k+6]*((datin.clm$TypeCity==1) + (datin.clm$TypeCounty==1) + (datin.clm$TypeTown==1) + (datin.clm$TypeVillage==1))
    z2 <- parms[k+7] + parms[k+8]*datin.clm$Coverage/1000 + 
          parms[k+9]*((datin.clm$TypeCity==1) + (datin.clm$TypeCounty==1) + (datin.clm$TypeTown==1) + (datin.clm$TypeVillage==1))
    rho1 <- z1
    rho2 <- -exp(z2)
    
  # frequency|deductible
  u1 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2),family=2,par=rho.both,par2=parms.both[1])
  u1a <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2-1),family=2,par=rho.both,par2=parms.both[1])
  # severity|deductible
  u2 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
  # density severity|deductible
  temp1 <- dGB2(yy3,mu,sigma,alpha1,alpha2)
  temp2 <- BiCopPDF(pGB2(yy1,mu1,sigma1,alpha11,alpha21),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
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



# The final model is cop1 = Frank, cop2 = rotated Gubmel 
ini.joint <- c(gam3,1,0,0,1,0,0)
loglikJoint(ini.joint,fam=c(5,24))  #31704.09
zop <- nlminb(ini.joint,loglikJoint,fam=c(5,24),control=list(eval.max=500))
print(zop)

l2.both <- zop$objective
l.both <- l1.both + l2.both

l1.coverage

loglik <- c(-l.entity,-l.coverage,-l.both)
AIC <- c(2*l.entity+2*38,2*l.coverage+2*38,2*l.both+2*41)
BIC <- c(2*l.entity+log(4152)*38,2*l.coverage+log(4152)*38,2*l.both+log(4152)*41)

# Table S.4: conditional copula models
cbind(loglik,AIC,BIC)

