# Application: ratemaking
# Need to run main_dep_int.cpp to obtain "mean_dep_int.txt"
# Need to run main_ind_int.cpp to obtain "mean_ind_int.txt"
# Need to run Prediction_Discrete.R to obtain "mean_discrete_int.txt"
# reproduce Table 4, Table 5, Figure 5, Figure 6, Table S.8


library(ggplot2)
library(mgcv)
library(numDeriv)
library(dglm)
library(ald)
library(quantreg)
library(rcompanion)
library(MASS)
library(cplm)


rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')

mean_ind_int <- read.table("mean_ind_int.txt")
mean_dep_int <- read.table("mean_dep_int.txt")
mean_discrete_int <- read.table("mean_discrete_int.txt")


# take year = 2006-2009; use 2010 as hold-out
selectyear <- 2010
datout.pol <- poldat[which(poldat$Year%in%selectyear),]

rsample <- cbind(datout.pol$PolicyNum, mean_ind_int, mean_dep_int,mean_discrete_int)
colnames(rsample) <- c("PolicyNum","IndInt","DepInt","DiscreteInt")
rsample <- merge(datout.pol,rsample,by="PolicyNum")
da <- data.frame(yobs=log(1+rsample$Claim),Pooled=1,Exogenous=log(rsample$IndInt), Endogenous=log(rsample$DepInt),Experience=log(rsample$Premium))
g<-gini(loss = "yobs", score=c("Pooled","Exogenous","Endogenous","Experience","yobs"), data = da)

# reproduce Table 4
gini <- g@gini[1:4,1:4]
sd <- g@sd[1:4,1:4]
print(gini);print(sd)


# reproduce Figure 5
plot(g@lorenz$Exogenous[,1],g@lorenz$Exogenous[,"Endogenous"],type="l",xlab="Premium(%)",ylab="Loss(%)",lty=2,main="Underwriting")
lines(g@lorenz$Exogenous[,1],g@lorenz$Exogenous[,"yobs"],type="l",lty=1)
abline(0,1)
legend("topleft",c("Underwriting Profit","Maximum Possible Profit"),lty=c(2,1))

plot(g@lorenz$Endogenous[,1],g@lorenz$Endogenous[,"Exogenous"],type="l",xlab="Premium(%)",ylab="Loss(%)",lty=2,main="Ratemaking")
lines(g@lorenz$Endogenous[,1],g@lorenz$Endogenous[,"Pooled"],type="l",lty=1)
abline(0,1)
legend("bottomright",c("Underwriting Loss","Maximum Possible Loss"),lty=c(2,1))



# produce gini index for entity types
# city
datcity <- subset(rsample,TypeCity==1)
da1 <- data.frame(yobs=log(1+datcity$Claim),Exogenous=log(datcity$IndInt), Endogenous=log(datcity$DepInt))
g1 <- gini(loss = "yobs", score=c("Exogenous","Endogenous"), data = da1)
# county
datcounty <- subset(rsample,TypeCounty==1)
da2 <- data.frame(yobs=log(1+datcounty$Claim),Exogenous=log(datcounty$IndInt), Endogenous=log(datcounty$DepInt))
g2 <- gini(loss = "yobs", score=c("Exogenous","Endogenous"), data = da2)
# school
datschool <- subset(rsample,TypeSchool==1)
da3 <- data.frame(yobs=log(1+datschool$Claim),Exogenous=log(datschool$IndInt), Endogenous=log(datschool$DepInt))
g3<-gini(loss = "yobs", score=c("Exogenous","Endogenous"), data = da3)
# Town
dattown <- subset(rsample,TypeTown==1)
da4 <- data.frame(yobs=log(1+dattown$Claim),Exogenous=log(dattown$IndInt), Endogenous=log(dattown$DepInt))
g4<-gini(loss = "yobs", score=c("Exogenous","Endogenous"), data = da4)
# village
datvillage <- subset(rsample,TypeVillage==1)
da5 <- data.frame(yobs=log(1+datvillage$Claim),Exogenous=log(datvillage$IndInt), Endogenous=log(datvillage$DepInt))
g5<-gini(loss = "yobs", score=c("Exogenous","Endogenous"), data = da5)
# misc
datmisc <- subset(rsample,TypeMisc==1)
da6 <- data.frame(yobs=log(1+datmisc$Claim),Exogenous=log(datmisc$IndInt), Endogenous=log(datmisc$DepInt))
g6<-gini(loss = "yobs", score=c("Exogenous","Endogenous"), data = da6)



# reproduce Figure 6
plot(g1@lorenz$Exogenous[,1],g1@lorenz$Exogenous[,"Endogenous"],type="l",xlab="Premium(%)",ylab="Loss(%)",main="City")
lines(g1@lorenz$Endogenous[,1],g1@lorenz$Endogenous[,"Exogenous"],type="l")
abline(0,1,lty=2)
text(x=30,y=80,"Underwriting Loss")
text(x=70,y=20,"Underwriting Profit")
arrows(30, 75, 50, 60)
arrows(70, 25, 65, 55)

plot(g2@lorenz$Exogenous[,1],g2@lorenz$Exogenous[,"Endogenous"],type="l",xlab="Premium(%)",ylab="Loss(%)",main="County")
lines(g2@lorenz$Endogenous[,1],g2@lorenz$Endogenous[,"Exogenous"],type="l")
abline(0,1,lty=2)
text(x=30,y=80,"Underwriting Loss")
text(x=70,y=20,"Underwriting Profit")
arrows(30, 75, 50, 60)
arrows(70, 25, 65, 55)

plot(g3@lorenz$Exogenous[,1],g3@lorenz$Exogenous[,"Endogenous"],type="l",xlab="Premium(%)",ylab="Loss(%)",main="School")
lines(g3@lorenz$Endogenous[,1],g3@lorenz$Endogenous[,"Exogenous"],type="l")
abline(0,1,lty=2)
text(x=30,y=80,"Underwriting Loss")
text(x=70,y=20,"Underwriting Profit")
arrows(30, 75, 50, 65)
arrows(70, 25, 65, 45)

plot(g4@lorenz$Exogenous[,1],g4@lorenz$Exogenous[,"Endogenous"],type="l",xlab="Premium(%)",ylab="Loss(%)",main="Town")
lines(g4@lorenz$Endogenous[,1],g4@lorenz$Endogenous[,"Exogenous"],type="l")
abline(0,1,lty=2)
text(x=30,y=80,"Underwriting Loss")
text(x=75,y=20,"Underwriting Profit")
arrows(30, 75, 35, 50)
arrows(75, 25, 55, 35)

plot(g5@lorenz$Exogenous[,1],g5@lorenz$Exogenous[,"Endogenous"],type="l",xlab="Premium(%)",ylab="Loss(%)",main="Village")
lines(g5@lorenz$Endogenous[,1],g5@lorenz$Endogenous[,"Exogenous"],type="l")
abline(0,1,lty=2)
text(x=30,y=80,"Underwriting Loss")
text(x=75,y=20,"Underwriting Profit")
arrows(30, 75, 35, 50)
arrows(75, 25, 50, 35)

plot(g6@lorenz$Exogenous[,1],g6@lorenz$Exogenous[,"Endogenous"],type="l",xlab="Premium(%)",ylab="Loss(%)",main="Miscellaneous")
lines(g6@lorenz$Endogenous[,1],g6@lorenz$Endogenous[,"Exogenous"],type="l")
abline(0,1,lty=2)
text(x=30,y=80,"Underwriting Loss")
text(x=75,y=20,"Underwriting Profit")
arrows(30, 75, 35, 58)
arrows(75, 25, 60, 38)

# reproduce Table 5
print(g1);print(g2);print(g3);print(g4);print(g5);print(g6)


# results in Table S.8: Discrete Deductible
da_disc <- data.frame(yobs=log(1+rsample$Claim),Exogenous=log(rsample$IndInt),Experience=log(rsample$Premium), DiscreteDeductible=log(rsample$DiscreteInt))
g_disc <- gini(loss = "yobs", score=c("Exogenous","Experience","DiscreteDeductible"), data = da_disc)
gini <- g_disc@gini[1:2,3]
sd <- g_disc@sd[1:2,3]
cbind(gini,sd)


