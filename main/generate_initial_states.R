# start 
rm(list = ls())
setwd("/Users/vincentarumadri/Desktop/Epi/Modelling/master_thesis_ua")

# packages 
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)

pacman::p_load(here, deSolve, tidyverse, patchwork, knitr, kableExtra, grid, gridExtra)

# From Hodgson et al 2020 
load(here::here("data/uk_data_sum.RData"))
load(here::here("data/posteriors.Rda"))
source("functions/RunInterventions.R")

# parameter values for analysis
param_means <- colMeans(post)

# using mean of @ parameter for analysis
currentParamValues <- param_means
ParameterValuesforODE(currentParamValues)

# generate inital M
ageStratification = uk_data_sum$ageGroupBoundary
populationPerAgeGroup = uk_data_sum$populationAgeGroup

initialStates_i <- generateInitialStates(param_means, ageStratification, populationPerAgeGroup)
initial_states <- matrix(as.numeric(initialStates_i), nrow = 26, byrow = FALSE)

# initial states list of each age group 
state_initial <- list()


for (i in 1:25) {
  current_initial_states <- as.numeric(initial_states[, i])
  
  state_initial[[i]] <- current_initial_states
  
}

state_initial
