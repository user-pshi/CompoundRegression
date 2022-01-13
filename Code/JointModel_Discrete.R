# joint model for claim frequency, claim severity, and discrete deductible choice
# estimate parameters in the severity model and copulas (R,Y) and (N|R,Y|R)
# Need to run DeductibleModels_Discrete.R to obtain "coef.deductClogit.RDS"
# Need to run FrequencyModels.R to obtain "coef.freq.RDS"
# Need to run Joint_FreqDeduct_Discrete.R to obtain "cop.DN_Discrete.RDS"
# reproduce results in Table S.7


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
# joint model of frequency and severity conditional on deductible
# use claim file
#########################################################################

datin.clm$deduct <- 1*(datin.clm$Deduct==500) + 2*(datin.clm$Deduct==1000) + 3*(datin.clm$Deduct==2500) +
  4*(datin.clm$Deduct==5000) + 5*(datin.clm$Deduct==10000) + 6*(datin.clm$Deduct==15000) +
  7*(datin.clm$Deduct==25000) + 8*(datin.clm$Deduct==50000) + 9*(datin.clm$Deduct==75000) +
  10*(datin.clm$Deduct==100000)
table(datin.clm$deduct)


# response variable
yy1 <- datin.clm$deduct
yy2 <- datin.clm$Freq
yy3 <- datin.clm$LossBeforeDeductible

# read in parameters
gam1 <- readRDS(file="coef.deductClogit.RDS")
gam2 <- readRDS(file="coef.freq.RDS")

copDN <- readRDS(file="cop.DN_Discrete.RDS")



######################
# deductible
######################
# covariates
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage)-1,data=datin.clm))

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
  #ans <- pmax(1e-10,ans)
  ans
}
# cdf of 01 inflated NB
Fcount <- function(y){
  ans <- pzero*(y>=0)+pone*(y>=1)+(1-pzero-pone)*pnbinom(y,mu=meannb,size=size)
  ans
}


######################
# severity
######################

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



# now define likelihood for the joint model

loglikJoint <- function(parms,fam){
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  rho1 <- parms[k+4]
  rho2 <- parms[k+5]
  df1 <- parms[k+6]
  df2 <- parms[k+7]
  
  # frequency|deductible
  u1 <- (BiCopCDF(Fclogit(yy1),Fcount(yy2),family=copDN[1],par=copDN[2],par2=copDN[3]) -
           BiCopCDF(Fclogit(yy1-1),Fcount(yy2),family=copDN[1],par=copDN[2],par2=copDN[3]))/fclogit(yy1)
  u1a <- (BiCopCDF(Fclogit(yy1),Fcount(yy2-1),family=copDN[1],par=copDN[2],par2=copDN[3]) -
            BiCopCDF(Fclogit(yy1-1),Fcount(yy2-1),family=copDN[1],par=copDN[2],par2=copDN[3]))/fclogit(yy1)
  # severity|deductible
  u2 <- (BiCopCDF(Fclogit(yy1),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1) - 
           BiCopCDF(Fclogit(yy1-1),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1))/fclogit(yy1)
  # density severity|deductible
  temp1 <- dGB2(yy3,mu,sigma,alpha1,alpha2)
  temp2 <- BiCopHfunc2(Fclogit(yy1),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1) -
    BiCopHfunc2(Fclogit(yy1-1),pGB2(yy3,mu,sigma,alpha1,alpha2),fam[1],rho1,df1)
  # (frequency,severity)|deductible
  u1 <- ifelse(u1<=0,1e-12,u1)
  u1a <- ifelse(u1a<=0,1e-12,u1a)
  u2 <- ifelse(u2<=0,1e-12,u2)
  u1 <- ifelse(u1>=1,1-1e-12,u1)
  u1a <- ifelse(u1a>=1,1-1e-12,u1a)
  u2 <- ifelse(u2>=1,1-1e-12,u2)
  temp3 <- BiCopHfunc2(u1,u2,fam[2],rho2,df2)-BiCopHfunc2(u1a,u2,fam[2],rho2,df2)
  
  # - log likelihood
  llk <- -sum(log(temp1)) - sum(log(temp2)) - sum(log(temp3))
  llk
}



# The final model is cop1 = student t, cop2 = rotated Gubmel 

ini.joint <- c(gam3,0.5,-1.5,5,5)
loglikJoint(ini.joint,fam=c(2,24))
zop <- nlminb(ini.joint,loglikJoint,fam=c(2,24),lower =c(rep(-Inf,ncol(cova)),-Inf,0.001,0.001,-50,-50),upper=c(rep(Inf,ncol(cova)),Inf,Inf,Inf,50,-1.001),control=list(eval.max=500))
print(zop)
npar <- length(zop$par)
hess<-hessian(loglikJoint,zop$par,fam=c(2,24))
se <-sqrt(diag(solve(hess[1:(npar-1),1:(npar-1)])))


# Estimated copula parameter (R,Y) and (N|R,Y|R) in Table S.7
# select fam = c(2,24)
cbind(zop$par[1:(npar-1)],se)[11:12,]

# Kenall's tau
# deduct and severity fam = 2 student t
BiCopPar2Tau(2,zop$par[11],zop$par[13])
# frequency and severity fam =24 Gumbel
BiCopPar2Tau(24,zop$par[12])


# save optimal parameter
saveRDS(c(fam=c(2,24),zop$par[1:13]),file="cop.DS_discrete.RDS")




