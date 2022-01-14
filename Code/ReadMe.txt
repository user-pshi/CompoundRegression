The codes can be used to perform analysis done in the paper titled "Copula Regression for Compound Distributions with Endogenous Covariates with Applications in Insurance Deductible Pricing" by Peng Shi and Gee Lee.

The programs are written in R and C++. Below are the table of contents for the codes:


I. R codes: descriptive statistics and estimation

SummaryStats.R: generate descriptive statistics (reproduce Figure 1, Figure 2, Figure 3, and Table 1)
DeductibleModels.R: fit models for the continuous deductible ratio (reproduce Table 2 and Figure S.1)
FrequencyModels.R: fit models for claim frequency (reproduce Table 2)
Joint_FreqDeduct.R: fit joint model for claim frequency and continuous deductible ratio (reproduce Table 2)
FrequencyComparison.R: Goodness-of-fit for marginal and conditional models for claim frequency (reproduce Table S.1)
JointModel.R: fit joint model for claim frequency, claim severity, and continuous deductible ratio (reproduce Table 2 and Figure S.1)
Prediction_FreqSev.R: perform hold-out sample prediction for claim frequency and claim severity (reproduce Table 3)
Application_a.R: perform analysis in the ratemaking application (reproduce Figure 5, Figure 6, Table 4, Table 5, Table S.8)
Application_b.R: perform analysis in the risk capital application (reproduce Table 6)

ModelComparison_a.R: model comparison for independence, conditional independence, and full dependence models (reproduce Table S.3 and Table S.4)
ModelComparison_b.R: model comparison for simplified and conditional copula models (reproduce Table S.4)

JointModel_Gaussian_a.R: analysis using Gaussian copula with Poisson frequency and Gamma severity (reproduce Table S.5) 
JointModel_Gaussian_b.R: analysis using Gaussian copula with Negative Binomial frequency and Lognormal severity (reproduce Table S.5) 
JointModel_Gaussian_c.R: analysis using Gaussian copula with ZOINB frequency and GB2 severity (reproduce Table S.5) 

JointModel_LossCost.R: fit joint model for continuous deductible and loss cost (reproduce Table S.6 and S.7)
Prediction_LossCost.R: perform hold-out sample prediction for loss cost (reproduce Table S.8)

DeductibleModels_Discrete.R: fit models for the discrete deductible choice (reproduce Table S.9)
Joint_FreqDeduct_Discrete.R: fit joint model for claim frequency and discrete deductible choice (reproduce Table S.7)
JointModel_Discrete.R: fit joint model for claim frequency, claim severity, and discrete deductible choice (reproduce Table S.7)
Prediction_Discrete.R: perform hold-out sample prediction using discrete deductible model (generate output "mean_discrete_int.txt")

Rds2Txt.R: convert data in Rds to Txt which are used by the C++ routines. 

II. C++ codes: prediction and simulation

Note: In general, the C++ codes can be compiled on a Unix-like system using commands such as: <Compiler Name> <File Name> -o <Output File Name> -std=c++11 -I<Include Directory> -L<Library Directory> -lgsl For example, to compile main_ind_int.cpp on a MacOS Big Sur console, we had to type: c++ main_ind_int.cpp -o main_ind_int -std=c++11 -I/usr/local/include -L/usr/local/lib -lgsl (Please suppress warning messages if your compiler gives minor warnings during the compilation.)

main_dep_int.cpp: predict aggregate claims using numerical integration approach for the dependence model (generate output "mean_dep_int.txt")
main_ind_int.cpp: predict aggregate claims using numerical integration approach for the independence model (generate output "mean_ind_int.txt")

The C++ codes below generate outputs that are required to generate workspace "Application.RData":
TestHDdep.cpp: simulate claims for the d=5000 (high deductible) case for the dependence model.
TestHDind.cpp: simulate claims for the d=5000 (high deductible) case for the independence model.
TestLDdep.cpp: simulate claims for the d=500 (low deductible) case for the dependence model.
TestLDind.cpp: simulate claims for the d=500 (low deductible) case for the independence model.

