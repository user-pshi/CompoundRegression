# joint model for claim frequency, claim severity, and deductible ratio
# estimate parameters in the severity model and copulas (R,Y) and (N|R,Y|R)
# Need to run DeductibleModels.R to obtain "coef.deductGB2.RDS"
# Need to run FrequencyModels.R to obtain "coef.freq.RDS"
# Need to run Joint_FreqDeduct.R to obtain "cop.DN.RDS"
# reproduce results in Table 2 and Figure S.1

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
# 
#########################################################################

datin.clm$dedratio <- datin.clm$Deduct/datin.clm$Coverage

# response variable
yy1 <- datin.clm$dedratio
yy2 <- datin.clm$Freq
yy3 <- datin.clm$LossBeforeDeductible

# read in parameters
gam1 <- readRDS(file="coef.deductGB2.RDS")
gam2 <- readRDS(file="coef.freq.RDS")
copDN <- readRDS(file="cop.DN.RDS")


# covariates for deductible
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

# save this estimate for independence model
saveRDS(zop$par,file="coef.sev.RDS")



# Define likelihood for the joint model

loglikJoint1 <- function(parms,fam){
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  rho1 <- parms[k+4]
  rho2 <- parms[k+5]
  df1 <- parms[k+6]
  df2 <- parms[k+7]
  
  # frequency|deductible
  u1 <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2),family=copDN[1],par=copDN[2],par2=copDN[3])
  u1a <- BiCopHfunc1(pGB2(yy1,mu1,sigma1,alpha11,alpha21),Fcount(yy2-1),family=copDN[1],par=copDN[2],par2=copDN[3])
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
ini.joint <- c(gam3,0.5,-1.5,5,5)
zop <- nlminb(ini.joint,loglikJoint,fam=c(5,24),lower =c(rep(-Inf,ncol(cova)),-Inf,0.001,0.001,-50,-50),upper=c(rep(Inf,ncol(cova)),Inf,Inf,Inf,50,-1.001),control=list(eval.max=500))
print(zop)
npar <- length(zop$par)
hess<-hessian(loglikJoint,zop$par,fam=c(5,24))
se <-sqrt(diag(solve(hess[1:(npar-2),1:(npar-2)])))

# reproduce Table 2: Severy component + copulas (R,Y) and (N|R,Y|R)
cbind(zop$par[1:(npar-2)],se)


# save optimal parameter
saveRDS(c(fam=c(5,24),zop$par[1:12]),file="cop.DS.RDS")



#############################################
# qq plot for severity model  Figure S.1
# qqplot based on conditional distribution
##############################################
qqnorm(loglikJoint1(zop$par,fam=c(5,24))$resid,xlim=c(-4,4),ylim=c(-4,4),cex=0.8,main="Claim Severity")
qqline(loglikJoint1(zop$par,fam=c(5,24))$resid)







################################
# Model selection
#################################
famset1 <- c(1,2,3,4,5,6,13,14,16)
famset2 <- c(1,2,5,23,24,26,33,34,36)
range1 <- matrix(c(-0.999,-0.999,0.001,1.001,-50,1.001,0.001,1.001,1.001,
                   0.999, 0.999,   50,   50, 50,   50,   50,   50,   50),length(famset1),2) 
range2 <- matrix(c(-0.999,-0.999,-50,   -50,   -50,   -50,   -50,   -50,   -50,
                   0.999, 0.999, 50,-0.001,-1.001,-1.001,-0.001,-1.001,-1.001),length(famset2),2) 
rho1 <- c(0.5,0.5,0.5,1.5,0.5,1.5,0.5,1.5,1.5)
rho2 <- c(-0.5,-0.5,-0.5,-0.5,-1.5,-1.5,-0.5,-1.5,-1.5)

sink("copulaselect.txt")
for (i in 1:length(famset1)){
  for (j in 1:length(famset2)){
    ini.joint <- c(gam3,rho1[i],rho2[j],5,5)
    zop <- nlminb(ini.joint,loglikJoint,fam=c(famset1[i],famset2[j]),lower =c(rep(-Inf,ncol(cova)),-Inf,0.001,0.001,range1[i,1],range2[j,1]),upper=c(rep(Inf,ncol(cova)),Inf,Inf,Inf,range1[i,2],range2[j,2]),control=list(eval.max=500))
    print(c(famset1[i],famset2[j]))
    print(zop)
  }
}
sink()


