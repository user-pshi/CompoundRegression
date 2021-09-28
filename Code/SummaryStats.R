rm(list=ls())
load(file='poldat.RData')
load(file='claimdat.RData')


# take year = 2006-2009; use 2010 as hold-out
selectyear <- 2006:2009
datin.pol <- poldat[which(poldat$Year%in%selectyear),]
datin.clm <- claimdat[which(claimdat$Year%in%selectyear),]


##################################################################################################
# plot of claim frequency verse individual severity
# use claim file
##################################################################################################

datin.clm2 <- datin.clm[-which(datin.clm$LossBeforeDeductible<=1),]
datin.clm2$FreqCat <- datin.clm2$Freq*(datin.clm2$Freq<=5) + 6*(datin.clm2$Freq>5)*(datin.clm2$Freq<=10)+7*(datin.clm2$Freq>10)
datin.clm2$FreqF <- factor(datin.clm2$FreqCat,labels=c("1","2","3","4","5","6-10",">10"))
datin.clm2$LossPerCov <- datin.clm2$LossBeforeDeductible/datin.clm2$Coverage


pdf("SevByFreq.pdf",height=6,width=9)
p <- ggplot(datin.clm2, aes(x=FreqF, y=LossPerCov)) + 
  geom_boxplot(outlier.size = 0, outlier.shape=NA, notch=TRUE) +
  #geom_violin()
  labs(title="Claim Severity by Frequency",x="Frequency", y = "Individual Severity") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10()  + theme_bw()                       
p
dev.off()


##################################################################################################
# distribution by entity types
############################################################################################
# severity
datin.clm2$EntityType <- 1*(datin.clm2$TypeCity==1) + 2*(datin.clm2$TypeCounty==1) + 3*(datin.clm2$TypeSchool==1) + 
  4*(datin.clm2$TypeTown==1) + 5*(datin.clm2$TypeVillage==1) + 6*(datin.clm2$TypeMisc==1)
datin.clm2$EntityType <- factor(datin.clm2$EntityType,labels=c("City","County","School","Town","Village","Others"))
table(datin.clm2$EntityType)

pdf("SevByEntity.pdf",height=6,width=6)
p <- ggplot(datin.clm2, aes(x=EntityType, y=LossPerCov)) + 
  geom_boxplot(outlier.size = 0, outlier.shape=NA, notch=TRUE) +
  labs(title="Claim Severity by Enity Type",x="Entity Type", y = "Individual Severity") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10()  + theme_bw()                      
p
dev.off()



##################################################################################################
# examine frequency and deductible
# use policy file
##################################################################################################

datin.pol2 <- datin.pol[-which(datin.pol$Freq==24),]
pdf("DeductByFreq.pdf",height=6,width=6)
p <- ggplot(datin.pol2, aes(x=Freq, y=dedratio)) + 
  geom_point() +
  labs(title="Deductible v.s. Claim Frequency",x="Frequency", y = "Deductible") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10() + theme_bw()  
p
dev.off()


##################################################################################################
# distribution by entity types
############################################################################################

datin.pol2$EntityType <- 1*(datin.pol2$TypeCity==1) + 2*(datin.pol2$TypeCounty==1) + 3*(datin.pol2$TypeSchool==1) + 
  4*(datin.pol2$TypeTown==1) + 5*(datin.pol2$TypeVillage==1) + 6*(datin.pol2$TypeMisc==1)
datin.pol2$EntityType <- factor(datin.pol2$EntityType,labels=c("City","County","School","Town","Village","Others"))
table(datin.pol2$EntityType)

# Deductible ratio
pdf("DeductByEntity.pdf",height=6,width=6)
p <- ggplot(datin.pol2, aes(x=EntityType, y=dedratio)) + 
  geom_boxplot(outlier.size = 0, outlier.shape=NA, notch=TRUE) +
  labs(title="Deductible by Entity Type",x="Entity Type", y = "Deductible") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10() + theme_bw()  
p
dev.off()

# claim frequency
pdf("FreqByEntity.pdf",height=6,width=7)
p <- ggplot(datin.pol2, aes(x=EntityType, y=Freq)) + 
  geom_count() +

  labs(title="Claim Frequency by Entity Type",x="Entity Type", y = "Frequency") + 
  theme(plot.title = element_text(hjust=0.5)) +
  theme_bw()  
p
dev.off()




##################################################################################################
# examine individual severity and deductible
# I need to use both policy and claim file, so that I get the zero claims from the policy file
##################################################################################################


datin.clm$dedratio <- datin.clm$Deduct/datin.clm$Coverage

clmdat1 <- data.frame(PolicyNum = datin.clm$PolicyNum, dedratio = datin.clm$dedratio, grounduploss = datin.clm$LossBeforeDeductible, Coverage = datin.clm$Coverage, gr=1)
clmdat2 <- data.frame(PolicyNum = datin.pol$PolicyNum[which(datin.pol$Freq==0)], dedratio = datin.pol$dedratio[which(datin.pol$Freq==0)], grounduploss = datin.pol$Claim[which(datin.pol$Freq==0)],Coverage = datin.pol$Coverage[which(datin.pol$Freq==0)],gr=0)
clmcombined <- rbind(clmdat1,clmdat2)  
clmcombined$LossPerCov <- clmcombined$grounduploss/clmcombined$Coverage
clmcombined$ClaimCount <- factor(clmcombined$gr,labels=c("=0",">0"))


pdf("DeductBySev.pdf",height=6,width=7)
p <- ggplot(clmcombined, aes(x=LossPerCov+1, y=dedratio),group=ClaimCount) + 
  geom_point(aes(shape=ClaimCount)) +
  #geom_violin()
  labs(title="Deductible v.s. Claim Severity",x="Individual Severity", y = "Deductible") + 
  theme(plot.title = element_text(hjust=0.5)) +
  scale_y_log10() + scale_x_log10() + theme_bw() 
p
dev.off()




