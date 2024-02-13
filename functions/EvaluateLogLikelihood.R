###
EvaluateLogLikelihood = function(dailyBirthRate_t, totPopulation_t, ageStratification_t){
  # initialization code
  A <<- length(ageStratification_t)
  populationPerAgeGroup <<- numeric(A)
  eta <<- c(0)
  for (i in 1:(A-1)) {
    populationPerAgeGroup[i] <<- dailyBirthRate_t*365*(ageStratification_t[i+1]-ageStratification_t[i])
    eta <<- c(eta, 1.0/(365.0*(ageStratification_t[i+1] - ageStratification_t[i])))
    modelIncidencePerTime <<- c(modelIncidencePerTime, 0)
  }
  modelIncidencePerTime <<- c(modelIncidencePerTime, 0)
  populationPerAgeGroup[A] <<- totPopulation_t - (dailyBirthRate_t*365)*ageStratification_t[A]
  eta <<- c(eta, dailyBirthRate_t/(totPopulation_t - (dailyBirthRate_t*365)*ageStratification_t[A]))
  
  # other initialization code
  dt <<- 1
  currentODETime <<- 0
  dayNoAfterBurn <<- 0
  weekNo <<- 0
  monthNo <<- 0
  valueLogLikelihood <<- 0
  loglikelihoodError <<- FALSE
}
  ## run_start = 0
  ## run_burn = 52*7 + 1
  ## run_oneyr = 52*7 + run_burn
  
  dt <- 1
  currentODETime <- 0
  dayNoAfterBurn <- 0
  weekNo <- 0
  valueLogLikelihood <- 0
  
  loglikelihoodError <- FALSE
  
###
  evaluateLogLikelihoodWeekly <- function(currentParamValues) {

  }
  
  evaluateLogLikelihoodMonthly <- function(currentParamValues) {
    
  }
  
  getWeeklySample <- function(currentParamValues, epFlag) {

  }
  
  getMonthlySampleCpp <- function(currentParamValues, epFlag) {
    # Implementation for getMonthlySampleCpp
  }
  
  getAnnualIncidenceCpp <- function(currentParamValues) {
    # Implementation for getAnnualIncidenceCpp
  }
  
  getProportionBornProtectedCpp <- function(currentParamValues) {
    # Implementation for getProportionBornProtectedCpp
  }
  
  # Define in R after construction
  contactMatrixPhy <- matrix(0, nrow = ..., ncol = ...)  # Replace ... with appropriate dimensions
  contactMatrixCon <- matrix(0, nrow = ..., ncol = ...)  # Replace ... with appropriate dimensions
  observedData <- matrix(0, nrow = ..., ncol = ...)  # Replace ... with appropriate dimensions
  
  lowerParamSupport <- numeric(length = 25)
  upperParamSupport <- numeric(length = 25)
  
  # Related to parameter values
  parameterValues <- numeric(25)
  
  parameterValuesforODE <- function(currentParamValues) {
    ep_t <- c(exp(currentParamValues[20] + seq(0, 15) * currentParamValues[21]), 
              rep(currentParamValues[22], times = 7), 
              rep(currentParamValues[23], times = 2))
    
    pA <- c(rep(currentParamValues[6], times = 12), 
            rep(currentParamValues[7], times = 4), 
            rep(currentParamValues[8], times = 2), 
            rep(currentParamValues[9], times = 7))
    
    list(ep_t = ep_t, pA = pA)
  }
  
