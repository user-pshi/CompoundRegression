# Application: risk capital
# reproduce Table 6


##################################################################################################
# Note:
# The output from C++ codes (TestHDdep.cpp, TestHDind.cpp, TestLDdep.cpp, TestLDind.cpp) are too large to store.
# We save the workspace "Application.RData" for replication purposes
# Readers should skip to line 71 to generate Table 6
# Interested readers should run C++ codes to generate outputs that are required to generate the workspace
#################################################################################################################

rm(list=ls())
load("Application.RData")


################################################################################################################################################################################
# This chunk require output from C++ codes: TestHDdep.cpp, TestHDind.cpp, TestLDdep.cpp, TestLDind.cpp
# This chunk will create the workspace "Application.RData"
# So skip this part to line 71 to generate table 6
################################################################################################################################################################################

# Helper functions
CI95 <- function(v) { c(mean(v)-qnorm(0.975)*sqrt(var(v)/length(v)), mean(v), mean(v)+qnorm(0.975)*sqrt(var(v)/length(v))) }
VAR <- function(v) { quantile(v,0.95) }
CTE <- function(v) { mean(v[v>quantile(v,0.95)]) }


# The following code requires the simulation output to be in directories simoutput_dep and simoutput_ind.
B <- 1000
Nrep <- 250
npol <- 1200
M_dep_all <- list()
M_ind_all <- list()
for (rep in 1:Nrep) {
  print(rep)
  M_dep <- matrix(NA,npol,B)
  M_ind <- matrix(NA,npol,B)
  for (i in 1:npol) {
    M_dep[i,] <- read.csv(paste("./simoutput_dep/",rep,"_",i,".txt",sep=""),header=FALSE)[,1]
    M_ind[i,] <- read.csv(paste("./simoutput_ind/",rep,"_",i,".txt",sep=""),header=FALSE)[,1]
  }
  M_dep_all <- cbind(M_dep_all,list(M_dep))
  M_ind_all <- cbind(M_ind_all,list(M_ind))
}

VARDepHD <- NULL; for(rep in 1:Nrep) { SampDepHD <- apply(M_dep_all[[rep]][  1: 600,],2,sum); VARDepHD <- c(VARDepHD,quantile(SampDepHD,0.95)) }
VARDepLD <- NULL; for(rep in 1:Nrep) { SampDepLD <- apply(M_dep_all[[rep]][601:1200,],2,sum); VARDepLD <- c(VARDepLD,quantile(SampDepLD,0.95)) }
VARDepPO <- NULL; for(rep in 1:Nrep) { SampDepPO <- apply(M_dep_all[[rep]]           ,2,sum); VARDepPO <- c(VARDepPO,quantile(SampDepPO,0.95)) }
VARIndHD <- NULL; for(rep in 1:Nrep) { SampIndHD <- apply(M_ind_all[[rep]][  1: 600,],2,sum); VARIndHD <- c(VARIndHD,quantile(SampIndHD,0.95)) }
VARIndLD <- NULL; for(rep in 1:Nrep) { SampIndLD <- apply(M_ind_all[[rep]][601:1200,],2,sum); VARIndLD <- c(VARIndLD,quantile(SampIndLD,0.95)) }
VARIndPO <- NULL; for(rep in 1:Nrep) { SampIndPO <- apply(M_ind_all[[rep]]           ,2,sum); VARIndPO <- c(VARIndPO,quantile(SampIndPO,0.95)) }
CTEDepHD <- NULL; for(rep in 1:Nrep) { SampDepHD <- apply(M_dep_all[[rep]][  1: 600,],2,sum); CTEDepHD <- c(CTEDepHD,mean(SampDepHD[SampDepHD>VARDepHD[rep]])) }
CTEDepLD <- NULL; for(rep in 1:Nrep) { SampDepLD <- apply(M_dep_all[[rep]][601:1200,],2,sum); CTEDepLD <- c(CTEDepLD,mean(SampDepLD[SampDepLD>VARDepLD[rep]])) }
CTEDepPO <- NULL; for(rep in 1:Nrep) { SampDepPO <- apply(M_dep_all[[rep]]           ,2,sum); CTEDepPO <- c(CTEDepPO,mean(SampDepPO[SampDepPO>VARDepPO[rep]])) }
CTEIndHD <- NULL; for(rep in 1:Nrep) { SampIndHD <- apply(M_ind_all[[rep]][  1: 600,],2,sum); CTEIndHD <- c(CTEIndHD,mean(SampIndHD[SampIndHD>VARIndHD[rep]])) }
CTEIndLD <- NULL; for(rep in 1:Nrep) { SampIndLD <- apply(M_ind_all[[rep]][601:1200,],2,sum); CTEIndLD <- c(CTEIndLD,mean(SampIndLD[SampIndLD>VARIndLD[rep]])) }
CTEIndPO <- NULL; for(rep in 1:Nrep) { SampIndPO <- apply(M_ind_all[[rep]]           ,2,sum); CTEIndPO <- c(CTEIndPO,mean(SampIndPO[SampIndPO>VARIndPO[rep]])) }
RMDepHD <- CTEDepHD/(CTEDepHD+CTEDepLD)*CTEDepPO
RMDepLD <- CTEDepLD/(CTEDepHD+CTEDepLD)*CTEDepPO
RMIndHD <- CTEIndHD/(CTEIndHD+CTEIndLD)*CTEIndPO
RMIndLD <- CTEIndLD/(CTEIndHD+CTEIndLD)*CTEIndPO
LDRatDep <- CTEDepLD/(CTEDepHD+CTEDepLD)
LDRatInd <- CTEIndLD/(CTEIndHD+CTEIndLD)

# Save the workspace.
# save.image("Application.RData")
#################################################################################
##################################################################################


# Generate Table 6 (in thousand dollars)
rm(list=ls())
load("Application.RData")

# risk capital
cte.h.ind <- CI95(RMIndHD)/1000
cte.l.ind <- CI95(RMIndLD)/1000
cte.h.dep <- CI95(RMDepHD)/1000
cte.l.dep <- CI95(RMDepLD)/1000
cte.ind <- CI95(CTEIndPO)/1000
cte.dep <- CI95(CTEDepPO)/1000

# high deductible 
cte.h.ind[c(1,3)]; cte.h.dep[c(1,3)];
# low deductible
cte.l.ind[c(1,3)]; cte.l.dep[c(1,3)]; 
# portfolio
cte.ind[c(1,3)]; cte.dep[c(1,3)];


# capital allocation
pct.h.ind <- CI95(1-LDRatInd)
pct.l.ind <- CI95(LDRatInd)
pct.h.dep <- CI95(1-LDRatDep)
pct.l.dep <- CI95(LDRatDep)
pct.ind <- CI95(1-2*LDRatInd)
pct.dep <- CI95(1-2*LDRatDep)

# high deductible 
pct.h.ind[c(1,3)]; pct.h.dep[c(1,3)]
# low deductible
pct.l.ind[c(1,3)]; pct.l.dep[c(1,3)]
# portfolio
pct.ind[c(1,3)]; pct.dep[c(1,3)]





