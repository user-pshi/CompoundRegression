# Joint Model for Deductible and Loss cost
# Continuous deductible using GB2
# Loss cost using Tweedie

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
library(tweedie)


rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009; use 2010 as hold-out
selectyear <- 2006:2009
datin.pol <- poldat[which(poldat$Year%in%selectyear),]
datin.clm <- claimdat[which(claimdat$Year%in%selectyear),]




###########################################################################
# joint model deductible and loss cost
# use policy file
#########################################################################

datin.pol$dedratio <- datin.pol$Deduct/datin.pol$Coverage

# response variable
yy1 <- datin.pol$dedratio
yy2 <- datin.pol$y


# deductible
# read in parameters
gam1 <- readRDS(file="coef.deductGB2.RDS")

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




# loss cost
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.pol))
k <- ncol(cova)                        # number of predictors


##########################################################################
###########################################################################
# This is to get the parameter for loss cost and save it
# no need to run in future
# Tweedie's compound Poisson model
fit.lc <- glm(y~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),
              family=tweedie(var.power=1.5,link.power=0), data=datin.pol)
summary(fit.lc)

loglikM<-function(parms){ 
  mu=exp(cova%*%parms[1:k])
  p=parms[k+1]
  phi=parms[k+2]
  dt <- dtweedie(yy2, p, mu, phi)
  llk <- -sum(log(dt))
  llk
}

init.tw<-c(coef(fit.lc),1.5,summary(fit.lc)$dispersion)
loglikM(init.tw)

zop <- nlminb(init.tw,loglikM, lower =c(rep(-Inf,7),1.001,0.001),upper =c(rep(Inf,7),1.999,Inf))
zop

hess<-hessian(loglikM,zop$par)
se <-sqrt(diag(solve(hess)))
tratio <- zop$par/se
cbind(zop$par,se,tratio)

saveRDS(zop$par,file="coef.tweedie.LC.RDS")
##########################################################################
###########################################################################


gam2 <- readRDS(file="coef.tweedie.LC.RDS")

mu.tw <- exp(cova%*%gam2[1:k])
p=gam2[k+1]
phi=gam2[k+2]



# now define likelihood for the joint model

loglikJ <- function(parms,fam){
  rho <- parms[1]
  df <- parms[2]
  u1 <- pGB2(yy1,mu,sigma,alpha1,alpha2)
  u2 <- ptweedie(yy2, p, mu.tw, phi)
  id <- which(yy2==0)
  
  temp1 <- BiCopHfunc1(u1[id],u2[id],fam,rho,df)
  temp2 <- BiCopPDF(u1[-id],u2[-id],fam,rho,df)
  
  temp3 <- dGB2(yy1,mu,sigma,alpha1,alpha2)
  temp4 <- dtweedie(yy2, p, mu.tw, phi)
    
  llk <- -sum(log(temp1)) -sum(log(temp2)) - sum(log(temp3)) -sum(log(temp4[-id]))
  llk
}

# The final selected copula is Frank

# fam = 5 frank 
init <- c(-0.5,5)
zop3 <- nlminb(init,loglikJ,fam=5,lower =c(-100,1),upper=c(100,Inf),control=list(eval.max=500))
print(zop3)
BiCopPar2Tau(5,-1.443643)
# par = -1.443643
# obj = 42924.58
# tau = -0.1537675






