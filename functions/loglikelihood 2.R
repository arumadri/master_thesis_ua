###
EvaluateLogLikelihood <- setRefClass("EvaluateLogLikelihood",
                                     fields = list(
                                       A = "numeric",
                                       populationPerAgeGroup = list(),
                                       eta = list(),
                                       modelIncidencePerTime = list(),
                                       pA = list(),
                                       ep_t = list(),
                                       dailyBirthRate = "numeric",
                                       totPopulation = "numeric",
                                       run_start = "numeric",
                                       run_burn = "numeric",
                                       run_oneyr = "numeric",
                                       run_full = "numeric",
                                       dt = "numeric",
                                       ageStratification = "numeric",
                                       dayNoAfterBurn = "numeric",
                                       weekNo = "numeric",
                                       monthNo = "numeric",
                                       valueLogLikelihood = "numeric",
                                       currentODETime = "numeric",
                                       loglikelihoodError = "logical",
                                       contactMatrixPhy = list(),
                                       contactMatrixCon = list(),
                                       observedData = list(),
                                       lowerParamSupport = list(),
                                       upperParamSupport = list(),
                                       parameterValues = list()
                                     ),
                                     methods = list(
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
                                       },
                                       evaluateLogLikelihoodCppWeekly = function(currentParamValues) {
                                         # implementation code
                                         ll <<- 0
                                         if (dayNoAfterBurn == 0){
                                           for (a in 1:A)
                                             x0 <<- c(x0, rep(0, 23))
                                         }
                                         if (dayNoAfterBurn %% 7 == 0 && dayNoAfterBurn > 0) {
                                           for (a in 1:A) {
                                             modelIncidencePerTime[a] <<- x0[22 + 23*(a-1)]
                                             inc_tot[a] <<- inc_tot[a] + modelIncidencePerTime[a]
                                           }
                                           ll <<- evaluateTimeStepLogLikelihood(weekNo)
                                           if (loglikelihoodError == TRUE)
                                             return(0)
                                           weekNo <<- weekNo + 1
                                         }
                                         checkStability(x0)
                                         if (loglikelihoodError == TRUE)
                                           return(0)
                                         getWeeklyLikelihood(x0, inc_tot)
                                         if (loglikelihoodError == TRUE)
                                           return(0)
                                         getWeeklyIncidence(x0, sampleWeeklyIncidence, TRUE)
                                         if (loglikelihoodError == TRUE)
                                           return(0)
                                         dayNoAfterBurn <<- dayNoAfterBurn + 1
                                         return(ll)
                                       },
                                       evaluateLogLikelihoodCppMonthly = function(currentParamValues) {
                                         # implementation code
                                         ll <<- 0
                                         if (dayNoAfterBurn == 0){
                                           for (a in 1:A)
                                             x0 <<- c(x0, rep(0, 23))
                                         }
                                         if (dayNoAfterBurn %% 30 == 0 && dayNoAfterBurn > 0) {
                                           for (a in 1:A) {
                                             modelIncidencePerTime[a] <<- x0[22 + 23*(a-1)]
                                             inc_tot[a] <<- inc_tot[a] + modelIncidencePerTime[a]
                                           }
                                           ll <<- evaluateTimeStepLogLikelihood(monthNo)
                                           if (loglikelihoodError == TRUE)
                                             return(0)
                                           monthNo <<- monthNo + 1
                                         }
                                         checkStability(x0)
                                         if (loglikelihoodError == TRUE)
                                           return(0)
                                         getMonthlyLikelihood(x0, inc_tot)
                                         if (loglikelihoodError == TRUE)
                                           return(0)
                                         getMonthlyIncidence(x0, sampleMonthlyIncidence, TRUE)
                                         if (loglikelihoodError == TRUE)
                                           return(0)
                                         dayNoAfterBurn <<- dayNoAfterBurn + 1
                                         return(ll)
                                       },
                                       getWeeklySampleCpp = function(currentParamValues, epFlag) {
                                         # implementation code
                                         for (a in 1:NoAgeG) {
                                           for (j in 1:9) {
                                             sampleWeeklyIncidence[weekNo, 9*(a-1) + j] <<- x0[455*(a-1) + 72*6 + 8 + j] * ep_t[(a-1)]
                                           }
                                         }
                                         dayNoAfterBurn <<- dayNoAfterBurn + 1
                                       },
                                       getMonthlySampleCpp = function(currentParamValues, epFlag) {
                                         # implementation code
                                         for (a in 1:A) {
                                           sampleMonthlyIncidence[monthNo, a] <<- x0[22 + 23*(a-1)] * ep_t[(a-1)]
                                         }
                                         dayNoAfterBurn <<- dayNoAfterBurn + 1
                                       },
                                       getAnnualIncidenceCpp = function(currentParamValues) {
                                         # implementation code
                                         if (dayNoAfterBurn == 0){
                                           for (a in 1:A)
                                             x0 <<- c(x0, rep(0, 23))
                                         }
                                         if (dayNoAfterBurn %% 7 == 0 && dayNoAfterBurn > 0) {
                                           for (a in 1:A) {
                                             annualIncidence[a] <<- annualIncidence[a] + x0[22 + 23*(a-1)]
                                           }
                                         }
                                         checkStability(x0)
                                         dayNoAfterBurn <<- dayNoAfterBurn + 1
                                       },
                                       getProportionBornProtectedCpp = function(currentParamValues) {
                                         # implementation code
                                         for (j in 19:21) {
                                           num_vec_wcba[j-18] <<- (x[1 + j*23] + x[6 + j*23] + x[11 + j*23] + x[16 + j*23] + x[2 + j*23] + x[7 + j*23] + x[12 + j*23] + x[17 + j*23])/populationPerAgeGroup[j]
                                         }
                                         sum_wcb <<- sum(num_vec_wcba)/3.0
                                         return(sum_wcb)
                                       },
                                       ParameterValuesforODE = function(currentParamValues) {
                                         # implementation code
                                         ep_t <<- list()
                                         pA <<- list()
                                         parameterValuesTemp <<- rep(0, 25)
                                         parameterValuesTemp["xi"] <<- currentParamValues[1]
                                         parameterValuesTemp["si"] <<- currentParamValues[2]
                                         parameterValuesTemp["ga0"] <<- currentParamValues[3]
                                         parameterValuesTemp["g1"] <<- currentParamValues[4]
                                         parameterValuesTemp["g2"] <<- currentParamValues[5]
                                         parameterValuesTemp["om"] <<- currentParamValues[6]
                                         parameterValuesTemp["pA1"] <<- currentParamValues[7]
                                         parameterValuesTemp["pA2"] <<- currentParamValues[8]
                                         parameterValuesTemp["pA3"] <<- currentParamValues[9]
                                         parameterValuesTemp["pA4"] <<- currentParamValues[10]
                                         parameterValuesTemp["alpha_i"] <<- currentParamValues[11]
                                         parameterValuesTemp["d1"] <<- currentParamValues[12]
                                         parameterValuesTemp["d2"] <<- currentParamValues[13]
                                         parameterValuesTemp["d3"] <<- currentParamValues[14]
                                         parameterValuesTemp["phi"] <<- currentParamValues[15]
                                         parameterValuesTemp["qp"] <<- currentParamValues[16]
                                         parameterValuesTemp["qc"] <<- currentParamValues[17]
                                         parameterValuesTemp["b1"] <<- currentParamValues[18]
                                         parameterValuesTemp["psi"] <<- currentParamValues[19]
                                         parameterValuesTemp["c5ep1"] <<- currentParamValues[20]
                                         parameterValuesTemp["c5ep2"] <<- currentParamValues[21]
                                         parameterValuesTemp["ep5"] <<- currentParamValues[22]
                                         parameterValuesTemp["ep6"] <<- currentParamValues[23]
                                         parameterValuesTemp["I1"] <<- currentParamValues[24]
                                         parameterValuesTemp["I2"] <<- currentParamValues[25]
                                         parameterValues <<- parameterValuesTemp
                                  
##
for (a in 1:16) {
  ep_t <<- c(ep_t, list(exp(parameterValuesTemp["c5ep1"] + (a-1)*parameterValuesTemp["c5ep2"])))
}
for (a in 17:22) {
  ep_t <<- c(ep_t, list(exp(parameterValuesTemp["ep5"]))))
}
for (a in 23:25) {
  ep_t <<- c(ep_t, list(exp(parameterValuesTemp["ep6"]))))
}
                                         
for (a in 1:12) {
  pA <<- c(pA, list(exp(parameterValuesTemp["pA1"]))))
}
for (a in 13:16) {
  pA <<- c(pA, list(exp(parameterValuesTemp["pA2"]))))
}
for (a in 17:18) {
  pA <<- c(pA, list(exp(parameterValuesTemp["pA3"]))))
}
for (a in 19:25) {
  pA <<- c(pA, list(exp(parameterValuesTemp["pA4"]))))
}

####
generateInitialStates = function() {
  # implementation code
  initial_M <<- c(rep(0, 607))
  for (a in 1:A) {
    x0 <<- c(x0, initial_M)
  }
},
initial_M = function() {
  # implementation code
  initialProportionExposure <<- c()
  for (a in 1:A) {
    initialProportionExposure <<- c(initialProportionExposure, ((parameterValues$alpha_i[a]*parameterValues$om)/(parameterValues$alpha_i[a]*parameterValues$om + parameterValues$g2[a])))
  }
  return(initialProportionExposure)
},
poisson_cdf = function(meann, y) {
  # implementation code
  cumMass <<- ppois(y, meann)
  return(cumMass)
},
evaluateTimeStepLogLikelihood = function(intervalNumber) {
  # implementation code
  ll_total <<- 0
  for (a in 1:A) {
    ll_exp <<- -ep_t[(a-1)]*x0[22 + 23*(a-1)] + modelIncidencePerTime[a] + log(ep_t[(a-1)]*x0[22 + 23*(a-1)])
    ll_total <<- ll_total + ll_exp
  }
  return(ll_total)
},
stirlingApproximation = function(n) {
  # implementation code
  approx_factorial <<- sqrt(2*pi*n)*(n/exp(1))^n
  return(approx_factorial)
},
checkStability = function(x_check) {
  # implementation code
  for (age in 1:A) {
    if (x_check[23*(age-1) + 1] < 0 || x_check[23*(age-1) + 2] < 0 || x_check[23*(age-1) + 3] < 0 || x_check[23*(age-1) + 4] < 0 || x_check[23*(age-1) + 5] < 0 || x_check[23*(age-1) + 6] < 0 || x_check[23*(age-1) + 7] < 0 || x_check[23*(age-1) + 8] < 0 || x_check[23*(age-1) + 9] < 0 || x_check[23*(age-1) + 10] < 0 || x_check[23*(age-1) + 11] < 0 || x_check[23*(age-1) + 12] < 0 || x_check[23*(age-1) + 13] < 0 || x_check[23*(age-1) + 14] < 0 || x_check[23*(age-1) + 15] < 0 || x_check[23*(age-1) + 16] < 0 || x_check[23*(age-1) + 17] < 0 || x_check[23*(age-1) + 18] < 0 || x_check[23*(age-1) + 19] < 0){
      loglikelihoodError <<- TRUE
      return()
    }
  }
},
check_incidence = function(incidence) {
  # implementation code
  for (a in 1:A) {
    if (is.na(incidence[a]) || incidence[a] < 0 || is.infinite(incidence[a])) {
      loglikelihoodError <<- TRUE
      return()
    }
  }
},
getWeeklyLikelihood = function(x_state, inc_tota) {
  # implementation code
  for (a in 1:A) {
    if (is.na(inc_tot[a])) {
      loglikelihoodError <<- TRUE
      return()
    }
  }
},
getWeeklyIncidence = function(x_state, weeklyIncidence, epFlag) {
  # implementation code
  for (a in 1:A) {
    if (is.na(sampleWeeklyIncidence[weekNo, 9*(a-1) + 1]) || sampleWeeklyIncidence[weekNo, 9*(a-1) + 1] < 0 || is.infinite(sampleWeeklyIncidence[weekNo, 9*(a-1) + 1])) {
      loglikelihoodError <<- TRUE
      return()
    }
  }
},
getMonthlyLikelihood = function(x_state, inc_total) {
  # implementation code
  for (a in 1:A) {
    if (is.na(inc_tot[a])) {
      loglikelihoodError <<- TRUE
      return()
    }
  }
},
getMonthlyIncidence = function(x_state, monthlyIncidence, epFlag) {
  # implementation code
  for (a in 1:A) {
    if (is.na(sampleMonthlyIncidence[monthNo, a]) || sampleMonthlyIncidence[monthNo, a] < 0 || is.infinite(sampleMonthlyIncidence[monthNo, a])) {
      loglikelihoodError <<- TRUE
      return()
    }
  }
}
)
)
###

generateInitialStates <- function() {
  initialStates <- c()
  
  populationMatPro <- initial_M()
  
  I1 <- parameterValues["I1"]
  I2 <- parameterValues["I2"]
  I3 <- 0.5
  
  si <- 1.0/parameterValues["si"]
  g0 <- 1.0/parameterValues["ga0"]
  g1 <- 1.0/((parameterValues["ga0"])*(parameterValues["g1"]))
  g2 <- 1.0/((parameterValues["ga0"])*(parameterValues["g1"])*(parameterValues["g2"]))
  d1 <- parameterValues["d1"]
  d2 <- parameterValues["d1"]*parameterValues["d2"]
  d3 <- parameterValues["d1"]*parameterValues["d2"]*parameterValues["d3"]
  
  for (a in 1:A) {
    if (a < A){
      a1 <- ageStratification(a)
      a2 <- ageStratification(a+1)
    } else {
      a1 <- ageStratification(a)
      a2 <- 90
    }
    
    propEachExposureGroup <- initialProportionExposure(I3, a1, a2)
    populationNoMatPro <- populationPerAgeGroup[a] - populationMatPro[a]
    
    initialStates <- c(initialStates, populationMatPro[a])
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[1]*(1-I1)*(1-I2))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[1]*I1*si/(si+g0))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[1]*I1*g0/(si+g0)*pA[a])
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[1]*I1*g0/(si+g0)*(1-pA[a]))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[1]*(1-I1)*I2)
    
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[2]*(1-d1*I1)*(1-I2))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[2]*d1*I1*si/(si+g1))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[2]*d1*I1*g1/(si+g1)*pA[a])
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[2]*d1*I1*g1/(si+g1)*(1-pA[a]))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[2]*(1-d1*I1)*I2)
    
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[3]*(1-d2*I1)*(1-I2))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[3]*d2*I1*si/(si+g2))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[3]*d2*I1*g2/(si+g2)*pA[a])
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[3]*d2*I1*g2/(si+g2)*(1-pA[a]))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[3]*(1-d2*I1)*I2)
    
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[4]*(1-d3*I1)*(1-I2))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[4]*d3*I1*si/(si+g2))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[4]*d3*I1*g2/(si+g2)*pA[a])
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[4]*d3*I1*g2/(si+g2)*(1-pA[a]))
    initialStates <- c(initialStates, populationNoMatPro*propEachExposureGroup[4]*(1-d3*I1)*I2)
    initialStates <- c(initialStates, 0)
    initialStates <- c(initialStates, 0)
  }
  
  return(initialStates)
}