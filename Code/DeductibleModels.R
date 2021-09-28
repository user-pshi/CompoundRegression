# Fit models for deductible ratio

library(ggplot2)
library(mgcv)
library(numDeriv)
library(dglm)
library(ald)
library(quantreg)
library(rcompanion)
library(MASS)


rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009; use 2010 as hold-out
selectyear <- 2006:2009
datin.pol <- poldat[which(poldat$Year%in%selectyear),]
datin.clm <- claimdat[which(claimdat$Year%in%selectyear),]




#########################################################################
# below I explored GB2 
#
######################################################################


# fit a gb2 distribution

# severity (GB2)
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

yy <- datin.pol$dedratio
cova <- model.matrix(fit1)

# define loglikelihood function
loglikM <- function(parms){         
  k <- ncol(cova)                        # number of predictors
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  llk <- -sum(log(dGB2(yy, mu, sigma, alpha1, alpha2)))
  resid <- qnorm(pGB2(yy, mu, sigma, alpha1, alpha2))
  ans <- list(llk=llk,resid=resid)
  ans
}
loglikM2 <- function(parms) loglikM(parms)$llk


init <- c(fit1$coefficients,1,1,1)
zop <- nlminb(init,loglikM2, lower =c(rep(-Inf,ncol(cova)),-Inf,0.001,0.001),control=list(eval.max=500,iter.max=500))   # make sure convergence
zop


saveRDS(zop$par,file="coef.deductGB2.RDS")



# loglik =  27376.89

# get standard error
hess<-hessian(loglikM2,zop$par)
se <-sqrt(diag(solve(hess)))
tratio <- zop$par/se
cbind(zop$par,se,tratio)


pdf("qq_deduct.pdf",width=6,height=6)
qqnorm(loglikM(zop$par)$resid,xlim=c(-4,4),ylim=c(-4,4),cex=0.8,main="Deductible")
qqline(loglikM(zop$par)$resid)
dev.off()


# histogram of residuals
gam <- readRDS("coef.deductGB2.RDS")
mu <- cova%*%gam[1:ncol(cova)]
sigma <- gam[ncol(cova)+1]
resid <- (log(yy)-mu)/sigma
pdf("hist_deductresid.pdf",width=6,height=6)
hist(resid,xlab="Deductible",main="")
dev.off()















