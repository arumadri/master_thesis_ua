# start 
rm(list = ls())
setwd("/Users/vincentarumadri/Desktop/Epi/Modelling/master_thesis_ua")
library(here)
# From Hodgson et al 2020 
load(here::here("data/uk_data_sum.RData"))
load(here::here("data/posteriors.Rda"))
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
initial_M(param_means, ageStratification, populationPerAgeGroup)

# # generate initial states
# pVHR <- uk_data_sum$pVHR
# pHR <- uk_data_sum$pHR
# pLR <- uk_data_sum$pLR
# p_mat <- uk_data_sum$prop_mat
# cov_c <- uk_data_sum$nmat
# initialStates <- generateInitialStates(cov_c, param_means, ageStratification, populationPerAgeGroup)
pVHR <- uk_data_sum$pVHR
pHR <- uk_data_sum$pHR
pLR <- uk_data_sum$pLR
p_mat <- uk_data_sum$prop_mat
cov_c <- uk_data_sum$nmat
initialStates_i <- generateInitialStates(p_mat, param_means, ageStratification, populationPerAgeGroup)

state_names <- c("M", "S0", "E0", "A0", "I0", "R0",
                      "S1", "E1", "A1", "I1", "R1",
                      "S2", "E2", "A2", "I2", "R2",
                      "S3", "E3", "A3", "I3", "R3"
                )
initialStates_matrix <- matrix(initialStates_i, nrow = 25, byrow = FALSE)

initialStates_df<- as.data.frame(initialStates_matrix)
colnames(initialStates_df) <- state_names

initialStates_df

### option 2

state_initial_values <- state_initial_age
state_initial_values[] = 0
# source [https://www.ethnicity-facts-figures.service.gov.uk/uk-population-by-ethnicity/demographics/age-groups/latest/]
state_initial_values[grepl('S0.', names(state_initial_values))] = c(3232055, 41846595, 14518785)
state_initial_values[grepl('I0.', names(state_initial_values))] = 1

params$age_r=0
time_values <- seq(0, (1*365), 1)
run_model = ode(y=state_initial_values, times = time_values, func = rsv_model_age, parms = params, method = "ode23")
plot(run_model[,'E01'])
plot(run_model[,'A01'])
plot(run_model[,'I01'])
plot(run_model[,'R01'])
plot(run_model[,'S02'])
head(run_model)
tail(run_model)
