# Descriptive Statistics
# reproduce results in Figure 1, Figure 2, Figure 3, and Table 1 


library(ggplot2)
library(ggpubr)

rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009 as training; use 2010 as hold-out
selectyear <- 2006:2009
datin.pol <- poldat[which(poldat$Year%in%selectyear),]
datin.clm <- claimdat[which(claimdat$Year%in%selectyear),]

# create additional variables in policy data for plotting
datin.pol$dedratio <- datin.pol$Deduct/datin.pol$Coverage
datin.pol$EntityType <- 1*(datin.pol$TypeCity==1) + 2*(datin.pol$TypeCounty==1) + 3*(datin.pol$TypeSchool==1) + 
  4*(datin.pol$TypeTown==1) + 5*(datin.pol$TypeVillage==1) + 6*(datin.pol$TypeMisc==1)
datin.pol$EntityType <- factor(datin.pol$EntityType,labels=c("City","County","School","Town","Village","Others"))

# create additional variables in claim data for plotting
datin.clm$dedratio <- datin.clm$Deduct/datin.clm$Coverage
datin.clm$EntityType <- 1*(datin.clm$TypeCity==1) + 2*(datin.clm$TypeCounty==1) + 3*(datin.clm$TypeSchool==1) + 
  4*(datin.clm$TypeTown==1) + 5*(datin.clm$TypeVillage==1) + 6*(datin.clm$TypeMisc==1)
datin.clm$EntityType <- factor(datin.clm$EntityType,labels=c("City","County","School","Town","Village","Others"))
datin.clm$FreqCat <- datin.clm$Freq*(datin.clm$Freq<=5) + 6*(datin.clm$Freq>5)*(datin.clm$Freq<=10)+7*(datin.clm$Freq>10)
datin.clm$FreqF <- factor(datin.clm$FreqCat,labels=c("1","2","3","4","5","6-10",">10"))
datin.clm$LossPerCov <- datin.clm$LossBeforeDeductible/datin.clm$Coverage




##################################################################################################
# Figure 1a: examine deductible and frequency 
##################################################################################################
p <- ggplot(datin.pol, aes(x=Freq, y=dedratio)) + 
  geom_point() +
  labs(title="Deductible v.s. Claim Frequency",x="Frequency", y = "Deductible") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10() + theme_bw()  
p


##################################################################################################
# Figure 1b: deductible and individual severity
# I need to use both policy and claim file, so that I get the zero claims from the policy file
##################################################################################################
clmdat1 <- data.frame(PolicyNum = datin.clm$PolicyNum, dedratio = datin.clm$dedratio, grounduploss = datin.clm$LossBeforeDeductible, Coverage = datin.clm$Coverage, gr=1)
clmdat2 <- data.frame(PolicyNum = datin.pol$PolicyNum[which(datin.pol$Freq==0)], dedratio = datin.pol$dedratio[which(datin.pol$Freq==0)], grounduploss = datin.pol$Claim[which(datin.pol$Freq==0)],Coverage = datin.pol$Coverage[which(datin.pol$Freq==0)],gr=0)
clmcombined <- rbind(clmdat1,clmdat2)  
clmcombined$LossPerCov <- clmcombined$grounduploss/clmcombined$Coverage
clmcombined$ClaimCount <- factor(clmcombined$gr,labels=c("=0",">0"))


p <- ggplot(clmcombined, aes(x=LossPerCov+1, y=dedratio),group=ClaimCount) + 
  geom_point(aes(shape=ClaimCount)) +
  #geom_violin()
  labs(title="Deductible v.s. Claim Severity",x="Individual Severity", y = "Deductible") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10() + scale_x_log10() + theme_bw() 
p


##################################################################################################
# Figure 2: claim frequency verse individual severity
##################################################################################################
p <- ggplot(datin.clm, aes(x=FreqF, y=LossPerCov)) + 
  geom_boxplot(outlier.size = 0, outlier.shape=NA, notch=TRUE) +
  #geom_violin()
  labs(title="Claim Severity by Frequency",x="Frequency", y = "Individual Severity") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10()  + theme_bw()                       
p




##################################################################################################
# Figure 3
############################################################################################
# Deductible by entity type
p <- ggplot(datin.pol, aes(x=EntityType, y=dedratio)) + 
  geom_boxplot(outlier.size = 0, outlier.shape=NA, notch=TRUE) +
  labs(title="Deductible by Entity Type",x="Entity Type", y = "Deductible") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10() + theme_bw()  
p


# claim frequency by entity type
p <- ggplot(datin.pol, aes(x=EntityType, y=Freq)) + 
  geom_count() +

  labs(title="Claim Frequency by Entity Type",x="Entity Type", y = "Frequency") + 
  theme(plot.title = element_text(hjust=0.5)) +
  theme_bw()  
p

# severity by entity type
p <- ggplot(datin.clm, aes(x=EntityType, y=LossPerCov)) + 
  geom_boxplot(outlier.size = 0, outlier.shape=NA, notch=TRUE) +
  labs(title="Claim Severity by Enity Type",x="Entity Type", y = "Individual Severity") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10()  + theme_bw()                      
p




##################################################################################################
# Generate Results in Table 1
############################################################################################
# response variables
mean(datin.pol$dedratio);sd(datin.pol$dedratio)
mean(datin.pol$Coverage);sd(datin.pol$Coverage)
mean(datin.pol$Freq);sd(datin.pol$Freq)

# rating variables
sum(datin.pol$TypeCity)/nrow(datin.pol)
sum(datin.pol$TypeCounty)/nrow(datin.pol)
sum(datin.pol$TypeSchool)/nrow(datin.pol)
sum(datin.pol$TypeTown)/nrow(datin.pol)
sum(datin.pol$TypeVillage)/nrow(datin.pol)
sum(datin.pol$TypeMisc)/nrow(datin.pol)
mean(datin.clm$LossBeforeDeductible);sd(datin.clm$LossBeforeDeductible)


