# Hold-out sample prediction using discrete deductible model
# Need to run DeductibleModels_Discrete.R to obtain "coef.deductClogit.RDS"
# Need to run FrequencyModels.R to obtain "coef.freq.RDS"
# Need to run Joint_FreqDeduct_Discrete.R to obtain "cop.DN_Discrete.RDS"
# Need to run JointModel_Discrete.R to obtain "cop.DS_discrete.RDS"
# Generates output mean_discrete_int.txt


library(ggplot2)
library(mgcv)
library(numDeriv)
library(ald)
library(quantreg)
library(MASS)
library(BB)
library(pscl)
library(VineCopula)
library(tweedie)
library(cplm)

rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')
selectyear <- 2006:2009
datin.pol <- poldat[which(poldat$Year%in%selectyear),]
datin.clm <- claimdat[which(claimdat$Year%in%selectyear),]
datin.clm$dedratio <- datin.clm$Deduct/datin.clm$Coverage

selectyearout <- 2010
datout.pol <- poldat[which(poldat$Year%in%selectyearout),]
datout.clm <- claimdat[which(claimdat$Year%in%selectyearout),]
datout.clm$dedratio <- datout.clm$Deduct/datout.clm$Coverage

# response variable
yy1 <- datout.pol$dedratio
yy2 <- datout.pol$Freq
yy3 <- datout.pol$LossBeforeDeductible

# read in parameters
gam1 <- readRDS(file="coef.deductCLogit.RDS")
gam2 <- readRDS(file="coef.freq.RDS")
copDN <- readRDS(file="cop.DN_discrete.RDS")
copDS <- readRDS(file="cop.DS_discrete.RDS")

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

# covariates
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datout.pol))
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


k <- 7
mu3dep <- cova%*%copDS[3:(k+2)]
sigma3dep <- copDS[k+3]
alpha13dep <- copDS[k+4]
alpha23dep <- copDS[k+5]
thetaNY <- copDS[2+k+3+2]
fam1 <- copDS[1]
fam2 <- copDS[2]

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

# pdf of 01 inflated NB
fcounti <- function(y,i){
  ans <- pzero[i]*(y==0)+pone[i]*(y==1)+(1-pzero[i]-pone[i])*dnbinom(y,mu=meannb[i],size=size)
  #ans <- pmax(1e-10,ans)
  ans
}
# cdf of 01 inflated NB
Fcounti <- function(y,i){
  ans <- pzero[i]*(y>=0)+pone[i]*(y>=1)+(1-pzero[i]-pone[i])*pnbinom(y,mu=meannb[i],size=size)
  ans
}

datout.pol$deduct <- 1*(datout.pol$Deduct==500) + 2*(datout.pol$Deduct==1000) + 3*(datout.pol$Deduct==2500) +
  4*(datout.pol$Deduct==5000) + 5*(datout.pol$Deduct==10000) + 6*(datout.pol$Deduct==15000) +
  7*(datout.pol$Deduct==25000) + 8*(datout.pol$Deduct==50000) + 9*(datout.pol$Deduct==75000) +
  10*(datout.pol$Deduct==100000)
table(datout.pol$deduct)

#######################################
# endogenous prediction
#######################################

yy <- datout.pol$deduct
pDeduct <- function(y) {
  parms <- gam1
  a <- parms[1:9]
  b <- parms[10:15]
  p1 <- exp(a[1] + cova[,-1]%*%b)/(1+exp(a[1] + cova[,-1]%*%b))
  p2 <- exp(a[2] + cova[,-1]%*%b)/(1+exp(a[2] + cova[,-1]%*%b))
  p3 <- exp(a[3] + cova[,-1]%*%b)/(1+exp(a[3] + cova[,-1]%*%b))
  p4 <- exp(a[4] + cova[,-1]%*%b)/(1+exp(a[4] + cova[,-1]%*%b))
  p5 <- exp(a[5] + cova[,-1]%*%b)/(1+exp(a[5] + cova[,-1]%*%b))
  p6 <- exp(a[6] + cova[,-1]%*%b)/(1+exp(a[6] + cova[,-1]%*%b))
  p7 <- exp(a[7] + cova[,-1]%*%b)/(1+exp(a[7] + cova[,-1]%*%b))
  p8 <- exp(a[8] + cova[,-1]%*%b)/(1+exp(a[8] + cova[,-1]%*%b))
  p9 <- exp(a[9] + cova[,-1]%*%b)/(1+exp(a[9] + cova[,-1]%*%b))
  pDeduct <- (y==1)*p1 + (y==2)*p2 + (y==3)*p3 + (y==4)*p4 + (y==5)*p5 + (y==6)*p6 +
    (y==7)*p7 + (y==8)*p8 + (y==9)*p9 + (y==10)*1
  pDeduct
}


dDeduct <- function(yy){
  parms <- gam1
  a <- parms[1:9]
  b <- parms[10:15]
  p1 <- exp(a[1] + cova[,-1]%*%b)/(1+exp(a[1] + cova[,-1]%*%b))
  p2 <- exp(a[2] + cova[,-1]%*%b)/(1+exp(a[2] + cova[,-1]%*%b))
  p3 <- exp(a[3] + cova[,-1]%*%b)/(1+exp(a[3] + cova[,-1]%*%b))
  p4 <- exp(a[4] + cova[,-1]%*%b)/(1+exp(a[4] + cova[,-1]%*%b))
  p5 <- exp(a[5] + cova[,-1]%*%b)/(1+exp(a[5] + cova[,-1]%*%b))
  p6 <- exp(a[6] + cova[,-1]%*%b)/(1+exp(a[6] + cova[,-1]%*%b))
  p7 <- exp(a[7] + cova[,-1]%*%b)/(1+exp(a[7] + cova[,-1]%*%b))
  p8 <- exp(a[8] + cova[,-1]%*%b)/(1+exp(a[8] + cova[,-1]%*%b))
  p9 <- exp(a[9] + cova[,-1]%*%b)/(1+exp(a[9] + cova[,-1]%*%b))
  dDeduct <- (yy==1)*p1 + (yy==2)*(p2-p1) + (yy==3)*(p3-p2) + (yy==4)*(p4-p3) + (yy==5)*(p5-p4) + (yy==6)*(p6-p5) +
    (yy==7)*(p7-p6)  + (yy==8)*(p8-p7) +  (yy==9)*(p9-p8) + (yy==10)*(1-p9)
  dDeduct
}


MaxN <- 20
nobs <- dim(cova)[1]
mean.dep<-rep(0,nobs)

for (i in 1:nobs){ 
  print(i)
  pp0 <- pzero[i]+(1-pzero[i]-pone[i])*dnbinom(0,mu=meannb[i],size=size)
  pp1 <- pzero[i]+(1-pzero[i]-pone[i])*dnbinom(0,mu=meannb[i],size=size)+pone[i]+(1-pzero[i]-pone[i])*dnbinom(1,mu=meannb[i],size=size)
  
  # Dependence model
  getmean_dep <- function(deduct) {
    
    y1 <- deduct 
    uHi <- pDeduct(y1) 
    uLo <- ifelse(pDeduct(y1-1)==0,0.00000001,pDeduct(y1-1))
    densNY_R_i <- function(y2,y3,i) { 
      uN_R  <- (BiCopCDF(uHi[i], Fcounti(y2,i),family=copDN[1],par=copDN[2],par2=copDN[3]) - BiCopCDF(uLo[i], Fcounti(y2,i),family=copDN[1],par=copDN[2],par2=copDN[3]))/dDeduct(y1)[i]
      uNa_R <- (BiCopCDF(uHi[i], ifelse(y2==0,0,Fcounti(y2-1,i)),family=copDN[1],par=copDN[2],par2=copDN[3]) - BiCopCDF(uLo[i], ifelse(y2==0,0,Fcounti(y2-1,i)),family=copDN[1],par=copDN[2],par2=copDN[3]))/dDeduct(y1)[i]
      if(uN_R>1) {
        temp <- uN_R
        uN_R <- uN_R * (1/temp/1.000000001)
        uNa_R <- uNa_R * (1/temp/1.000000001)
      }
      uY_R       <- (BiCopCDF(uHi[i], pGB2(y3,mu3dep[i],sigma3dep,alpha13dep,alpha23dep),family=copDS[1],par=copDS[13],par2=copDS[15]) -
                       BiCopCDF(uLo[i], pGB2(y3,mu3dep[i],sigma3dep,alpha13dep,alpha23dep),family=copDS[1],par=copDS[13],par2=copDS[15]))/dDeduct(y1)[i]
      uY_R_upper <- (BiCopCDF(uHi[i], pGB2(datout.pol[i,]$Coverage*1000000,mu3dep[i],sigma3dep,alpha13dep,alpha23dep),family=copDS[1],par=copDS[13],par2=copDS[15]) - 
                       BiCopCDF(uLo[i], pGB2(datout.pol[i,]$Coverage*1000000,mu3dep[i],sigma3dep,alpha13dep,alpha23dep),family=copDS[1],par=copDS[13],par2=copDS[15]))/dDeduct(y1)[i]
      if(uY_R_upper>1) {
        temp <- uY_R_upper
        uY_R_upper <- uY_R_upper * (1/temp/1.000000001)
        uY_R <- uY_R * (1/temp/1.000000001)
      }
      cdf_NY_R_Lo_i <- BiCopCDF(Fcounti(y2-1,i),uY_R_upper,fam2,thetaNY,1)
      cdf_NY_R_Hi_i <- BiCopCDF(Fcounti(y2,i),uY_R_upper,fam2,thetaNY,1)
      dens_Y_R_i <- (BiCopHfunc2(uHi[i],pGB2(y3,mu3dep[i],sigma3dep,alpha13dep,alpha23dep),family=copDS[1],par=copDS[13],par2=copDS[15]) -
                       BiCopHfunc2(uLo[i],pGB2(y3,mu3dep[i],sigma3dep,alpha13dep,alpha23dep),family=copDS[1],par=copDS[13],par2=copDS[15])) *
        dGB2(y3,mu3dep[i],sigma3dep,alpha13dep,alpha23dep) / dDeduct(y1)[i]
      densNY_R_i <- (BiCopHfunc2(uN_R,uY_R,fam2,thetaNY,1) - BiCopHfunc2(uNa_R,uY_R,fam2,thetaNY,1)) * dens_Y_R_i
      densNY_R_i
    }
    
    funcY <- function(y) { sum(sapply(1:MaxN,function(n){n*densNY_R_i(n,y,i)})) } 
    funcYvec <- function(yvec) {sapply(yvec,function(y){y*funcY(y)})}
    term1 <- integrate(funcYvec,lower=0,upper=datout.pol[i,]$Coverage*1000000)$value
    term2 <- datout.pol[i,]$Coverage*1000000*sum(sapply(1:MaxN,function(n){
      
      y1 <- deduct # discrete deductible
      uHi <- pDeduct(y1) 
      uLo <- ifelse(pDeduct(y1-1)==0,0.00000001,pDeduct(y1-1))
      y2 <- n
      uN_R <- (BiCopCDF(uHi[i], Fcounti(y2,i),family=copDN[1],par=copDN[2],par2=copDN[3]) - BiCopCDF(uLo[i], Fcounti(y2,i),family=copDN[1],par=copDN[2],par2=copDN[3]))/dDeduct(y1)[i]
      uNa_R <- (BiCopCDF(uHi[i], ifelse(y2==0,0,Fcounti(y2-1,i)),family=copDN[1],par=copDN[2],par2=copDN[3]) - BiCopCDF(uLo[i], ifelse(y2==0,0,Fcounti(y2-1,i)),family=copDN[1],par=copDN[2],par2=copDN[3]))/dDeduct(y1)[i]
      if(uN_R>1) {
        temp <- uN_R
        uN_R <- uN_R * (1/temp/1.000000001)
        uNa_R <- uNa_R * (1/temp/1.000000001)
      }
      upp <- datout.pol[i,]$Coverage*1000000
      uY_R <- (BiCopCDF(uHi[i], pGB2(upp,mu3dep[i],sigma3dep,alpha13dep,alpha23dep),family=copDS[1],par=copDS[13],par2=copDS[15])-
                 BiCopCDF(uLo[i], pGB2(upp,mu3dep[i],sigma3dep,alpha13dep,alpha23dep),family=copDS[1],par=copDS[13],par2=copDS[15]))/dDeduct(y1)[i]
      u1 <- BiCopCDF(uN_R,uY_R,fam2,thetaNY,1)
      u1a <- BiCopCDF(uNa_R,uY_R,fam2,thetaNY,1)
      densN_R_i <- (BiCopCDF(uHi[i],Fcounti(y2,i),family=copDN[1],par=copDN[2],par2=copDN[3]) - 
                      BiCopCDF(uLo[i],Fcounti(y2,i),family=copDN[1],par=copDN[2],par2=copDN[3]) - 
                      BiCopCDF(uHi[i],ifelse(y2==0,0,Fcounti(y2-1,i)),family=copDN[1],par=copDN[2],par2=copDN[3]) + 
                      BiCopCDF(uLo[i],ifelse(y2==0,0,Fcounti(y2-1,i)),family=copDN[1],par=copDN[2],par2=copDN[3])) / dDeduct(y1)[i]
      n * ( densN_R_i - (u1 - u1a) )
    }))
  }
  
  mean.dep[i] <- getmean_dep(datout.pol$deduct[i])
}


write.table(mean.dep, "mean_discrete_int.txt", col.names=FALSE, row.names=FALSE)


