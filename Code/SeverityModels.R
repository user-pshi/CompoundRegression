# Fit models for claim severity

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


###########################################################################
# Analyze severity: this is to perform a marginal analysis, the only purpose is to get initial estimates
# use claim file
#########################################################################

yy <- datin.clm$LossBeforeDeductible

cova <- model.matrix(lm(yy~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datin.clm))
k <- ncol(cova)  



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

# define loglikelihood function
loglikM <- function(parms){         
  mu <- cova%*%parms[1:k]
  sigma <- parms[k+1]
  alpha1 <- parms[k+2]
  alpha2 <- parms[k+3]
  llk <- -sum(log(dGB2(yy, mu, sigma, alpha1, alpha2)))
  llk
}



# save parameters
saveRDS(zop$par,file="coef.sev.RDS")

