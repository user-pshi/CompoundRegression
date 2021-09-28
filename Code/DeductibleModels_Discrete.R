# Fit models for discrete deductible

library(ggplot2)
library(mgcv)
library(numDeriv)
library(dglm)
library(ald)
library(quantreg)
library(rcompanion)
library(MASS)
library(VGAM)


rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009; use 2010 as hold-out
selectyear <- 2006:2009
datin.pol <- poldat[which(poldat$Year%in%selectyear),]
datin.clm <- claimdat[which(claimdat$Year%in%selectyear),]


###########################################################################
# Analyze deductible 
# use policy file
# use discrete choice model
#########################################################################
names(datin.pol)
table(datin.pol$Deduct)
datin.pol$deduct <- 1*(datin.pol$Deduct==500) + 2*(datin.pol$Deduct==1000) + 3*(datin.pol$Deduct==2500) +
                    4*(datin.pol$Deduct==5000) + 5*(datin.pol$Deduct==10000) + 6*(datin.pol$Deduct==15000) +
                    7*(datin.pol$Deduct==25000) + 8*(datin.pol$Deduct==50000) + 9*(datin.pol$Deduct==75000) +
                   10*(datin.pol$Deduct==100000)
table(datin.pol$deduct)

library(ordinal)
fit <- clm(factor(deduct)~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage), data=datin.pol)
summary(fit)


yy <- datin.pol$deduct
cova <- model.matrix(lm(yy~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage)-1,data=datin.pol))


# code likelihood
loglikD <- function(parms){
  a <- parms[1:9]
  b <- parms[10:15]
  p1 <- exp(a[1] + cova%*%b)/(1+exp(a[1] + cova%*%b))
  p2 <- exp(a[2] + cova%*%b)/(1+exp(a[2] + cova%*%b))
  p3 <- exp(a[3] + cova%*%b)/(1+exp(a[3] + cova%*%b))
  p4 <- exp(a[4] + cova%*%b)/(1+exp(a[4] + cova%*%b))
  p5 <- exp(a[5] + cova%*%b)/(1+exp(a[5] + cova%*%b))
  p6 <- exp(a[6] + cova%*%b)/(1+exp(a[6] + cova%*%b))
  p7 <- exp(a[7] + cova%*%b)/(1+exp(a[7] + cova%*%b))
  p8 <- exp(a[8] + cova%*%b)/(1+exp(a[8] + cova%*%b))
  p9 <- exp(a[9] + cova%*%b)/(1+exp(a[9] + cova%*%b))
  
  l <- (yy==1)*p1 + (yy==2)*(p2-p1) + (yy==3)*(p3-p2) + (yy==4)*(p4-p3) + (yy==5)*(p5-p4) + (yy==6)*(p6-p5) +
       (yy==7)*(p7-p6)  + (yy==8)*(p8-p7) +  (yy==9)*(p9-p8) + (yy==10)*(1-p9)
  llk <- -sum(log(l))
  llk
}

init <- c(coef(fit)[1:9],-coef(fit)[10:15])
loglikD(init)

zop <- nlminb(init,loglikD)
zop

hess<-hessian(loglikD,zop$par,method.args=list(d=0.01))
se <-sqrt(diag(solve(hess)))
tratio <- zop$par/se
cbind(zop$par,se,tratio)



saveRDS(zop$par,file="coef.deductCLogit.RDS")






