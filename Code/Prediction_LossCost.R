# Hold-out sample prediction for loss cost
# Need to run DeductibleModels.R to obtain "coef.deductGB2.RDS"
# Need to run JointModel_LossCost.R to obtain "cop.DLC.RDS"
# reproduce results in Table S.8

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
library(cplm)


rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009; use 2010 as hold-out
datin.pol <- poldat[which(poldat$Year%in%c(2006:2009)),]
datout.pol <- poldat[which(poldat$Year%in%c(2010)),]

# read in parameters in deductible
gam1 <- readRDS(file="coef.deductGB2.RDS")
copDLC <- readRDS(file="cop.DLC.RDS")

# estimate parameters in Tweedie model
fit.lc <- glm(Claim~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),
              family=tweedie(var.power=1.5,link.power=0), data=datin.pol)
cova <- model.matrix(lm(Claim~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.pol))
k <- ncol(cova) 
loglikM<-function(parms){ 
  mu=exp(cova%*%parms[1:k])
  p=parms[k+1]
  phi=parms[k+2]
  dt <- dtweedie(datin.pol$Claim, p, mu, phi)
  llk <- -sum(log(dt))
  llk
}

init.tw<-c(coef(fit.lc),1.5,summary(fit.lc)$dispersion)
loglikM(init.tw)
zop <- nlminb(init.tw,loglikM, lower =c(rep(-Inf,7),1.001,0.001),upper =c(rep(Inf,7),1.999,Inf))
gam2 <- zop$par

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


###########################################################################
# Calculate predicted value
#########################################################################

# response variable
datout.pol$dedratio <- datout.pol$Deduct/datout.pol$Coverage
yy1 <- datout.pol$dedratio
yy2 <- datout.pol$Claim


# deductible
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage,data=datout.pol))
k1 <- ncol(cova1)

mu <- cova1%*%gam1[1:k1]
sigma <- gam1[k1+1]
alpha1 <- gam1[k1+2]
alpha2 <- gam1[k1+3]


# loss cost
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datout.pol))
k <- ncol(cova)                        # number of predictors

mu.tw <- exp(cova%*%gam2[1:k])
p=gam2[k+1]
phi=gam2[k+2]


# exogenous prediction
lc.ex <- mu.tw

# endogenous prediction
h1 <- function(u,v){
  par = copDLC[2]
  temp1 <- exp(-par*u)*(exp(-par*v)-1)
  temp2 <- exp(-par)-1 + (exp(-par*u)-1)*(exp(-par*v)-1)
  ans <- temp1/temp2
  ans
}

intg <- function(s,k){
  u1 <- pGB2(yy1[k],mu[k],sigma,alpha1,alpha2)
  u2 <- ptweedie(s, p, mu.tw[k], phi)
  ans <- 1-h1(u1,u2)
  ans
}

predf <- function(k) integrate(intg,lower=0,upper=5e7,k=k)$value
vpredf <- Vectorize(predf)
k <- 1:length(yy2)
lc.en <- vpredf(k)

da <- data.frame(yobs=log(1+yy2),Exogenous=log(1+lc.ex), Experience=log(datout.pol$Premium), Endogenous=log(1+lc.en))
g<-gini(loss = "yobs", score=c("Exogenous","Experience","Endogenous"), data = da)

# results in Table S.8: Aggregate Loss
gini <- g@gini[1:2,3]
sd <- g@sd[1:2,3]
cbind(gini,sd)
