Insurance claims data used in the paper titled "Copula Regression for Compound Distributions with Endogenous Covariates with Applications in Insurance Deductible Pricing" by Peng Shi and Gee Lee.


poldat.RData: policy-level information for the sample of policyholders
claimdat.RData: claim-level information for the sample of policyholders


Variable Description:

PolicyNum: policy ID
ClaimNum: claim ID
Year: policy year
Deduct: per-occurrence deductible
Coverage: amount of coverage
Premium: annual premium
Freq: number of claims
Claim: total amount of claims
LossBeforeDeductible: group-up losses
EntityType: 
	TypeCity: indicator of city 
	TypeCounty: indicator of county
	TypeSchool: indicator of school district 
	TypeTown: indicator of town 
	TypeVillage: indicator of village 
	TypeMisc: indicator of others


In addition, this folder contains some data which are results of intermediate steps. Some part of the analysis is time-consuming due to computational complexity. As a result, we save these intermediate results so that each code file could be run independently. Interested readers could skip certain parts to reproduce results in the paper.

mean_ind_int.txt: prediction from the exogenous deductible model 
mean_dep_int.txt: prediction from the endogenous deductible model (deductible ratio)
mean_discrete_int: prediction from the endogenous deductible model (deductible level)
Application.RData: simulated loss cost in risk captial application

Other files in the folder are intermediate results that are used as input for the C++ codes.