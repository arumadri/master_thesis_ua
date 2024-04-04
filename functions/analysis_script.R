# start 
rm(list = ls())
setwd("/Users/vincentarumadri/Desktop/Epi/Modelling/master_thesis_ua")
source("functions/RunInterventions.R")

# parameter values for analysis
param_means <- colMeans(post)
param_medians <- apply(post, 2, median)

# using mean of @ parameter for analysis
currentParamValues <- param_means
ParameterValuesforODE(currentParamValues)

# generate inital M
ageStratification = uk_data_sum$ageGroupBoundary
populationPerAgeGroup = uk_data_sum$populationAgeGroup
initial_M(parameterValues, ageStratification, populationPerAgeGroup)

# generate initial states
pVHR <- uk_data_sum$pVHR
pHR <- uk_data_sum$pHR
pLR <- uk_data_sum$pLR
p_mat <- uk_data_sum$prop_mat
cov_c <- uk_data_sum$nmat
initialStates <- generateInitialStates(cov_c)

# debug error 
initialProportionExposure(0.5,75,90)
