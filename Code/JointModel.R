# joint model for claim frequency, claim severity, and deductible ratio

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

datin.clm$dedratio <- datin.clm$Deduct/datin.clm$Coverage

# response variable
yy1 <- datin.clm$dedratio
yy2 <- datin.clm$Freq
yy3 <- datin.clm$LossBeforeDeductible

# read in parameters
gam1 <- readRDS(file="coef.deductGB2.RDS")
gam2 <- readRDS(file="coef.freq.RDS")
gam3 <- readRDS(file="coef.sev.RDS")

copDN <- readRDS(file="cop.DN.RDS")

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
  #ans <- pmax(1e-10,ans)
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



# now define likelihood for the joint model

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




# try different copulas
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




# The final model is cop1 = Frank, cop2 = rotated Gubmel 
# re-estimate the model

ini.joint <- c(gam3,0.5,-1.5,5,5)
zop <- nlminb(ini.joint,loglikJoint,fam=c(5,24),lower =c(rep(-Inf,ncol(cova)),-Inf,0.001,0.001,-50,-50),upper=c(rep(Inf,ncol(cova)),Inf,Inf,Inf,50,-1.001),control=list(eval.max=500))
print(zop)
npar <- length(zop$par)

# par
#(Intercept)      TypeCity    TypeCounty    TypeSchool      TypeTown   TypeVillage log(Coverage)                                           
#7.52355887   -0.24870376   -0.53656562    0.24082766   -0.14434654   -0.11259569    0.02311667    0.86280015    2.09233744    0.52640255 
#1.07003085   -1.50412594    5.00000000    5.00000000 
# loglik = 31709.09


# save optimal parameter
saveRDS(c(fam=c(5,24),zop$par[1:12]),file="cop.DS.RDS")



