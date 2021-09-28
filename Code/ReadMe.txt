The codes can be used to perform analysis done in the paper titled "Copula Regression for Compound Distributions with Endogenous Covariates with Applications in Insurance Deductible Pricing" by Peng Shi and Gee Lee.

The programs are written in R and C++. Below are the table of contents for the codes:


I. R codes: descriptive statistics and estimation

SummaryStats.R: generate descriptive statistics 
FrequencyModels.R: fit models for claim frequency
SeverityModels.R: fit models for claim severity
DeductibleModels.R: fit models for the continuous deductible ratio
Joint_FreqDeduct.R: fit joint model for claim frequency and deductible
Joint_SevDeduct.R: fit joint model for claim severity and deductible
JointModel.R: fit joint model for claim frequency, claim severity, and deductible
JointModel_Gaussian.R: fit joint model for claim frequency, claim severity, and deductible using Gaussian copulas
JointModel_LossCost.R: fit joint model for loss cost (Tweedie) and deductible


II. C++ codes: prediction and simulation

main_cind_int2.cpp: Predict aggregatre claims using the integral approach for the conditional independence model.
main_cind_sim_freq.cpp: Predict frequencies using the simulation approach for the conditional independence model.
main_cind_sim.cpp: Predict aggregate claims using the simulation approach for the conditional independence model.
main_dep_int2.cpp: Predict aggregate claims using the integral approach for the dependence model.
main_dep_sim_freq.cpp: Predict frequencies using the simulation approach for the dependence model.
main_dep_sim.cpp: Predict aggregate claims using the simulation approach for the dependence model.
main_ind_int2.cpp: Predict aggregate claims using the integral approach for the independence model.
main_ind_sim_freq.cpp: Predict frequencies using the simulation approach for the independence model.
main_ind_sim.cpp: Predict aggregate claims using the simulation approach for the independence model.
main_gaussian_int2.cpp: Predict aggregate claims using the integral approach for the Gaussian model.
TestHDdep.cpp: Simulate claims for the d=5000 (high deductible) case for the dependence model.
TestHDind.cpp: Simulate claims for the d=5000 (high deductible) case for the independence model.
TestLDdep.cpp: Simulate claims for the d=500 (low deductible) case for the dependence model.
TestLDind.cpp: Simulate claims for the d=500 (low deductible) case for the independence model.
TestHDdep.cpp: Simulate claims for the d=5000 (high deductible) case for the Gaussian model.
TestLDdep.cpp: Simulate claims for the d=500 (low deductible) case for the Gaussian model.

