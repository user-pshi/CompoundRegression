The codes can be used to perform analysis done in the paper titled "Copula Regression for Compound Distributions with Endogenous Covariates with Applications in Insurance Deductible Pricing" by Peng Shi and Gee Lee.

The programs are written in R and C++. Below are the table of contents for the codes:


I. R codes: descriptive statistics and estimation

SummaryStats.R: generate descriptive statistics 
DeductibleModels.R: fit models for the continuous deductible ratio 
FrequencyModels.R: fit models for claim frequency 
Joint_FreqDeduct.R: fit joint model for claim frequency and continuous deductible ratio 
JointModel.R: fit joint model for claim frequency, claim severity, and continuous deductible ratio 
Prediction_FreqSev.R: perform hold-out sample prediction for claim frequency and claim severity 
Application_a.R: perform analysis in the ratemaking application 
Application_b.R: perform analysis in the risk capital application 
ModelComparison_a.R: model comparison for independence, conditional independence, and full dependence models 
ModelComparison_b.R: model comparison for simplified and conditional copula models 
Rds2Txt.R: convert data in Rds to Txt which are used by the C++ routines. 

II. C++ codes: prediction and simulation

Note: In general, the C++ codes can be compiled on a Unix-like system using commands such as: <Compiler Name> <File Name> -o <Output File Name> -std=c++11 -I<Include Directory> -L<Library Directory> -lgsl For example, to compile main_ind_int.cpp on a MacOS Big Sur console, we had to type: c++ main_ind_int.cpp -o main_ind_int -std=c++11 -I/usr/local/include -L/usr/local/lib -lgsl (Please suppress warning messages if your compiler gives minor warnings during the compilation.)

main_dep_int.cpp: predict aggregate claims using numerical integration approach for the dependence model 
main_ind_int.cpp: predict aggregate claims using numerical integration approach for the independence model 

TestHDdep.cpp: simulate claims for the d=5000 (high deductible) case for the dependence model.
TestHDind.cpp: simulate claims for the d=5000 (high deductible) case for the independence model.
TestLDdep.cpp: simulate claims for the d=500 (low deductible) case for the dependence model.
TestLDind.cpp: simulate claims for the d=500 (low deductible) case for the independence model.

