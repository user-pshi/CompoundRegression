# rds2txt
# The code will convert rds data into txt format
# The txt files are input to the C++ routines
# Run this before runnning C++ routines

rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')
selectyearout <- 2010
datout.pol <- poldat[which(poldat$Year%in%selectyearout),]
datout.clm <- claimdat[which(claimdat$Year%in%selectyearout),]

y1 <- datout.pol$Deduct / datout.pol$Coverage
write.table(y1,"dat_y1.txt", col.names=FALSE, row.names=FALSE)
write.table(datout.pol$Coverage,"dat_coverage.txt", col.names=FALSE, row.names=FALSE)
write.table(datout.pol$Deduct,"dat_deduct.txt", col.names=FALSE, row.names=FALSE)

# response variable
datout.pol$dedratio <- datout.pol$Deduct/datout.pol$Coverage
yy1 <- datout.pol$dedratio
yy2 <- datout.pol$Freq
yy3 <- datout.pol$LossBeforeDeductible

# covariates for deductible
cova1 <- model.matrix(lm(yy1~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage,data=datout.pol))
write.table(cova1,"dat_cova1.txt", col.names=FALSE, row.names=FALSE)

# covariates for frequency and severity
cova <- model.matrix(lm(yy2~TypeCity+TypeCounty+TypeSchool+TypeTown+TypeVillage+log(Coverage),data=datout.pol))
write.table(cova,"dat_cova.txt", col.names=FALSE, row.names=FALSE)

# PolicyNum
write.table(datout.pol$PolicyNum,"dat_policynum.txt", col.names=FALSE, row.names=FALSE)

# New coefficients.
coef.deductCLogit <- readRDS("coef.deductCLogit.RDS")
write.table(coef.deductCLogit,"coef.deductCLogit.txt", col.names=FALSE, row.names=FALSE)
coef.deductGB2 <- readRDS("coef.deductGB2.RDS")
write.table(coef.deductGB2,"coef.deductGB2.txt", col.names=FALSE, row.names=FALSE)
coef.freq.NB <- readRDS("coef.freq.NB.RDS")
write.table(coef.freq.NB,"coef.freq.NB.txt", col.names=FALSE, row.names=FALSE)
coef.freq.Poisson <- readRDS("coef.freq.Poisson.RDS")
write.table(coef.freq.Poisson,"coef.freq.Poisson.txt", col.names=FALSE, row.names=FALSE)
coef.freq <- readRDS("coef.freq.RDS")
write.table(coef.freq,"coef.freq.txt", col.names=FALSE, row.names=FALSE)
coef.sev <- readRDS("coef.sev.RDS")
write.table(coef.sev,"coef.sev.txt", col.names=FALSE, row.names=FALSE)
cop.DLC <- readRDS("cop.DLC.RDS")
write.table(cop.DLC,"cop.DLC.txt", col.names=FALSE, row.names=FALSE)
cop.DN_Discrete <- readRDS("cop.DN_Discrete.RDS")
write.table(cop.DN_Discrete,"cop.DN_Discrete.txt", col.names=FALSE, row.names=FALSE)
cop.DN <- readRDS("cop.DN.RDS")
write.table(cop.DN,"cop.DN.txt", col.names=FALSE, row.names=FALSE)
cop.DS_discrete <- readRDS("cop.DS_discrete.RDS")
write.table(cop.DS_discrete,"cop.DS_discrete.txt", col.names=FALSE, row.names=FALSE)
cop.DS <- readRDS("cop.DS.RDS")
write.table(cop.DS,"cop.DS.txt", col.names=FALSE, row.names=FALSE)

