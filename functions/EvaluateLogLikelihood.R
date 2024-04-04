###
# Initialization
A <- 0 
populationPerAgeGroup <- numeric()
eta <-  numeric()
modelIncidencePerTime <- numeric()
pA <-  c()
ep_t <-  c()

dailyBirthRate <- 0 
totPopulation <- 0  
run_start <- 0       
run_burn <- 0        
run_oneyr <- 0       
run_full <- 0        
dt <- 0              

ageStratification <- numeric()  

dayNoAfterBurn <- 0
weekNo <- 0
monthNo <- 0
valueLogLikelihood <- 0

currentODETime <- 0
loglikelihoodError <- FALSE
###
EvaluateLogLikelihood <- function(dailyBirthRate_t, totPopulation_t, ageStratification_t) {
  A <- length(ageStratification_t)
  eta <- c(0)
  populationPerAgeGroup <- c()
  modelIncidencePerTime <- rep(0, A)  # Pre-allocate space for A elements
  pA <- c()  # will fill this later 
  ep_t <- c()  # will fill this later 
  
  for (i in 1:(A-1)) {
    populationPerAgeGroup <- c(populationPerAgeGroup, dailyBirthRate_t*365*(ageStratification_t[i+1] - ageStratification_t[i]))
    eta <- c(eta, 1.0/(365.0*(ageStratification_t[i+1] - ageStratification_t[i])))
  }
  
  # total population and daily birth rate
  populationPerAgeGroup <- c(populationPerAgeGroup, totPopulation_t - (dailyBirthRate_t*365)*ageStratification_t[A])
  eta <- c(eta, dailyBirthRate_t/(totPopulation_t - (dailyBirthRate_t*365)*ageStratification_t[A]))
  
  # Other 
  #run_start <-  0
  #run_burn <-  52*7 + 1
  #run_oneyr <-  52*7 + run_burn
  dt <- 1
  currentODETime <- 0
  dayNoAfterBurn <- 0
  weekNo <- 0
  monthNo <- 0
  valueLogLikelihood <- 0
  loglikelihoodError <- FALSE
  
  # return
  return(list(
    A = A,
    populationPerAgeGroup = populationPerAgeGroup,
    eta = eta,
    modelIncidencePerTime = modelIncidencePerTime,
    pA = pA,
    ep_t = ep_t,
    dailyBirthRate = dailyBirthRate_t,
    totPopulation = totPopulation_t,
    ageStratification = ageStratification_t,
    dt = dt,
    currentODETime = currentODETime,
    dayNoAfterBurn = dayNoAfterBurn,
    weekNo = weekNo,
    monthNo = monthNo,
    valueLogLikelihood = valueLogLikelihood,
    loglikelihoodError = loglikelihoodError
  ))
}
###
evaluateLogLikelihoodCppWeekly <- function(currentParamValues) {
  
}

evaluateLogLikelihoodCppMonthly <- function(currentParamValues) {
  # Implement calculation logic here
}

getWeeklySampleCpp <- function(currentParamValues, epFlag) {
  # Implement calculation logic here
  return(sampleWeekly) # Assume sampleWeekly is the resulting matrix
}

getMonthlySampleCpp <- function(currentParamValues, epFlag) {
  # Implement calculation logic here
  return(sampleMonthly) # Assume sampleMonthly is the resulting matrix
}

getAnnualIncidenceCpp <- function(currentParamValues) {
  # Implement calculation logic here
  return(annualIncidence) # Assume annualIncidence is the resulting vector
}

getProportionBornProtectedCpp <- function(currentParamValues) {
  # Implement calculation logic here
  return(proportionProtected) # Assume proportionProtected is the resulting vector
}

# Variable initializations
contactMatrixPhy <- matrix(nrow = 0, ncol = 0) # Adjust dimensions as needed
contactMatrixCon <- matrix(nrow = 0, ncol = 0) # Adjust dimensions as needed
observedData <- matrix(nrow = 0, ncol = 0) # Adjust dimensions as needed
lowerParamSupport <- numeric() # Initialize as empty, to be filled as needed
upperParamSupport <- numeric() # Initialize as empty, to be filled as needed
parameterValues <- numeric() # Placeholder, to be filled with actual parameter values
###
ParameterValuesforODE <- function(currentParamValues) {
  ep_t <- list()
  pA <- list()
  parameterValuesTemp <- rep(0, 25)
  names(parameterValuesTemp) <- c("xi", "si", "ga0", "g1", "g2", "om", "pA1", "pA2",
                                  "pA3", "pA4", "alpha_i", "d1", "d2", "d3", "phi",
                                  "qp", "qc", "b1", "psi", "c5ep1", "c5ep2", "ep5",
                                  "ep6", "I1", "I2")
  parameterValuesTemp <- as.list(currentParamValues)
  
  for (a in 1:16) {
    ep_t <- c(ep_t, exp(parameterValuesTemp$c5ep1 + (a-1)*parameterValuesTemp$c5ep2))
  }
  for (a in 17:22) {
    ep_t <- c(ep_t, parameterValuesTemp$ep5)
  }
  for (a in 23:25) {
    ep_t <- c(ep_t, parameterValuesTemp$ep6)
  }
  
  for (a in 1:12) {
    pA <- c(pA, parameterValuesTemp$pA1)
  }
  for (a in 13:16) {
    pA <- c(pA, parameterValuesTemp$pA2)
  }
  for (a in 17:18) {
    pA <- c(pA, parameterValuesTemp$pA3)
  }
  for (a in 19:25) {
    pA <- c(pA, parameterValuesTemp$pA4)
  }
  
  # return
  return(list(ep_t = ep_t, pA = pA, parameterValues = parameterValuesTemp))
}
###
generateInitialStates <- function() {
  initialStates <- c()
  
  # source necessary functions
  source("functions/RunInterventions.R")
  populationMatPro <- initial_M()
  
  I1 <- parameterValues["I1"]
  I2 <- parameterValues["I2"]
  I3 <- 0.5
  
  si <- 1.0 / parameterValues["si"]
  g0 <- 1.0 / parameterValues["ga0"]
  g1 <- 1.0 / (parameterValues["ga0"] * parameterValues["g1"])
  g2 <- 1.0 / (parameterValues["ga0"] * parameterValues["g1"] * parameterValues["g2"])
  d1 <- parameterValues["d1"]
  d2 <- parameterValues["d1"] * parameterValues["d2"]
  d3 <- parameterValues["d1"] * parameterValues["d2"] * parameterValues["d3"]
  
  for (a in 1:A) {
    if (a < A) {
      a1 <- ageStratification[a]
      a2 <- ageStratification[a + 1]
    } else {
      a1 <- ageStratification[a]
      a2 <- 90
    }
    
    propEachExposureGroup <- initialProportionExposure(I3, a1, a2)
    populationNoMatPro <- populationPerAgeGroup[a] - populationMatPro[a]
    
    initialStates <- c(initialStates,
                       populationMatPro[a],
                       populationNoMatPro * propEachExposureGroup[1] * (1 - I1) * (1 - I2),
                       populationNoMatPro * propEachExposureGroup[1] * I1 * si / (si + g0),
                       populationNoMatPro * propEachExposureGroup[1] * I1 * g0 / (si + g0) * pA[a],
                       populationNoMatPro * propEachExposureGroup[1] * I1 * g0 / (si + g0) * (1 - pA[a]),
                       populationNoMatPro * propEachExposureGroup[1] * (1 - I1) * I2,
                       
                       populationNoMatPro * propEachExposureGroup[2]*(1-d1*I1)*(1-I2),
                       populationNoMatPro * propEachExposureGroup[2]*d1*I1*si/(si+g1),
                       populationNoMatPro * propEachExposureGroup[2]*d1*I1*g1/(si+g1)*pA[a],
                       populationNoMatPro * propEachExposureGroup[2]*d1*I1*g1/(si+g1)*(1-pA[a]),
                       populationNoMatPro * propEachExposureGroup[2]*(1-d1*I1)*I2,
                       
                       populationNoMatPro * propEachExposureGroup[3]*(1-d2*I1)*(1-I2),
                       populationNoMatPro * propEachExposureGroup[3]*d2*I1*si/(si+g2),
                       populationNoMatPro * propEachExposureGroup[3]*d2*I1*g2/(si+g2)*pA[a],
                       populationNoMatPro * propEachExposureGroup[3]*d2*I1*g2/(si+g2)*(1-pA[a]),
                       populationNoMatPro  *propEachExposureGroup[3]*(1-d2*I1)*I2,
                       
                       populationNoMatPro * propEachExposureGroup[4]*(1-d3*I1)*(1-I2),
                       populationNoMatPro * propEachExposureGroup[4]*d3*I1*si/(si+g2),
                       populationNoMatPro * propEachExposureGroup[4]*d3*I1*g2/(si+g2)*pA[a],
                       populationNoMatPro  *propEachExposureGroup[4]*d3*I1*g2/(si+g2)*(1-pA[a]),
                       populationNoMatPro * propEachExposureGroup[4]*(1-d3*I1)*I2,
                       0,
                       0)
   
  }
  
  return(initialStates)
}
###
initial_M <- function(xi, ageStratification, populationPerAgeGroup) {
  A <- length(ageStratification)
  # initialize storage vector
  init_con <- numeric(A) 
  
  for (i in 1:(A - 1)) {
    # proportion of the population in each age group
    init_con_temp <- (pexp(365 * ageStratification[i + 1], rate = xi) - pexp(365 * ageStratification[i], rate = xi)) / ((365 * ageStratification[i + 1] - 365 * ageStratification[i]) * xi)
    init_con[i] <- init_con_temp * populationPerAgeGroup[i]
  }
  
  # last age group separately
  init_con[A] <- (pexp(365 * 90, rate = xi) - pexp(365 * ageStratification[A], rate = xi)) / ((365 * 90 - 365 * ageStratification[A]) * xi) * populationPerAgeGroup[A]
  
  return(init_con)
}
###
initialProportionExposure <- function(l, a1, a2) {
  # Initialize a vector for proportions
  prop <- numeric(A)
  
    prop <- vector("numeric", 4)
    prop[1] <- abs(ppois(0, a2*l) - ppois(0, a1*l))/((a2-a1)*l)
    prop[2] <- abs(ppois(1, a2*l) - ppois(1, a1*l))/((a2-a1)*l)
    prop[3] <- abs(ppois(2, a2*l) - ppois(2, a1*l))/((a2-a1)*l)
    prop[4] <- 1 - (prop[3] + prop[2] + prop[1])
  
  return(prop)
}
### 
poisson_cdf <- function(l, a, x) {
  if (l == 0.0 || a == 0.0){
    p <- dpois(0:1, lambda = 0.000001)
    return(ppois(x, lambda = 0.000001))
  }
  else {
    p <- dpois(0:1, lambda = l*a)
    return(ppois(x, lambda = l*a))
  }
}
###
evaluateTimeStepLogLikelihood <- function(metricNo, A, observedData, modelIncidencePerTime, ep_t) {
  ll <- 0
  loglikelihoodError <- FALSE
  
  for (a in 1:A) {
    dataNewInfections <- observedData[metricNo, a]
    if (dataNewInfections > modelIncidencePerTime[a]) {
      loglikelihoodError <- TRUE
      return(log(0))
    }
    
    estimatedLogBinomialCoeff <- stirlingApproximation(modelIncidencePerTime[a]) - 
      stirlingApproximation(dataNewInfections) - 
      stirlingApproximation(modelIncidencePerTime[a] - dataNewInfections)
    ll <- ll + estimatedLogBinomialCoeff + 
      dataNewInfections * log(ep_t[a]) + 
      (modelIncidencePerTime[a] - dataNewInfections) * log(1 - ep_t[a])
  }
  
  if (is.infinite(ll) || is.nan(ll)) {
    loglikelihoodError <- TRUE
  }
  
  return(list(logLikelihood = ll, logLikelihoodError = loglikelihoodError))
}
### 
stirlingApproximation <- function(n) {
  if (n == 0) {
    return(0)
  } else {
    x <- n + 1
    return((x - 0.5) * log(x) - x + 0.5 * log(2 * pi) + 1.0 / (12 * x) - 1.0 / (360 * x^3))
  }
}
###
checkStability <- function(x, A) {
  loglikelihoodError <- FALSE
  X_w <- 0
  
  # iteration limit
  iterationLimit <- A * 23
  
  if(length(x) != iterationLimit) {
    stop("The length of vector x does not match the expected structure based on A.")
  }
  
  # if any state values are negative
  if(any(x[1:iterationLimit] < 0)) {
    loglikelihoodError <- TRUE
    return(list(loglikelihoodError = loglikelihoodError, message = "Negative values found."))
  }
  
  X_w <- sum(x[1:iterationLimit])
  
  # if the sum of state values is non-finite (inf or NaN)
  if(is.infinite(X_w) || is.nan(X_w)) {
    loglikelihoodError <- TRUE
    return(list(loglikelihoodError = loglikelihoodError, message = "Non-finite sum detected."))
  }
  
  # loglikelihoodError status, should be FALSE
  return(list(loglikelihoodError = loglikelihoodError, message = "Stability check passed."))
}

###
getProportionBornProtected <- function(x, populationPerAgeGroup) {
  num_vec_wcba <- numeric(0) # Initialize an empty numeric vector
  sum_wcb <- 0.0
  
  for(j in 18:20) { 
    CB2_temp <- (x[2 + j*23] + x[7 + j*23] + x[12 + j*23] + x[17 + j*23] + 
                   x[3 + j*23] + x[8 + j*23] + x[13 + j*23] + x[18 + j*23]) / populationPerAgeGroup[j]
    num_vec_wcba <- c(num_vec_wcba, CB2_temp / 3.0)
  }
  
  sum_wcb <- sum(num_vec_wcba)
  
  return(sum_wcb)
}
###
check_incidence <- function(inc_tot, populationPerAgeGroup) {
  pl <- 0
  for (a in 1:25) {
    if (inc_tot[a] > populationPerAgeGroup[a] * 0.8) {
      return(log(0)) # Returns -Inf 
    }
  }
  return(pl) # Returns 0 if none of the conditions are met
}
###
getWeeklyLikelihood <- function(x0, inc_tot, dayNoAfterBurn, weekNo, A, loglikelihoodError, valueLogLikelihood) {
  if (dayNoAfterBurn == 0) {
    x0[seq(22, by = 23, length.out = A) + 1] <- 0 # Reset incidence at t_d = 0; Adjust for R's indexing
  }
  
  if (dayNoAfterBurn %% 7 == 0 && dayNoAfterBurn > 0) {
    for (a in 1:A) {
      index <- 22 + a * 23
      modelIncidencePerTime <- x0[index] # Extract model incidence for the week
      inc_tot[a] <- inc_tot[a] + modelIncidencePerTime
      x0[index] <- 0 # Reset for next period
    }
    # Update log likelihood for the week
    evaluateResult <- evaluateTimeStepLogLikelihood(weekNo) # Placeholder for actual function call
    valueLogLikelihood <- valueLogLikelihood + evaluateResult$valueLogLikelihood
    loglikelihoodError <- evaluateResult$loglikelihoodError || loglikelihoodError
    
    if (loglikelihoodError) {
      return(list(x0 = x0, inc_tot = inc_tot, dayNoAfterBurn = dayNoAfterBurn, weekNo = weekNo, loglikelihoodError = loglikelihoodError, valueLogLikelihood = valueLogLikelihood))
    }
    weekNo <- weekNo + 1
  }
  
  # Annual reset of inc_tot
  if (dayNoAfterBurn %% 365 == 0 && dayNoAfterBurn > 0) {
    # ll += check_incidence(inc_tot)
    inc_tot <- numeric(A) # Reset inc_tot for the new year
  }
  
  dayNoAfterBurn <- dayNoAfterBurn + 1
  
  return(list(x0 = x0, inc_tot = inc_tot, dayNoAfterBurn = dayNoAfterBurn, weekNo = weekNo, loglikelihoodError = loglikelihoodError, valueLogLikelihood = valueLogLikelihood))
}
###
getWeeklyIncidence <- function(x0, sampleWeeklyIncidence, epFlag, dayNoAfterBurn, weekNo, A, ep_t) {
  if (dayNoAfterBurn == 0) {
    x0[seq(23, by = 23, length.out = A) + 22] <- 0 # Reset incidence at t_d = 0
  }
  
  if (dayNoAfterBurn %% 7 == 0 && dayNoAfterBurn > 0) {
    for (a in 1:A) {
      for (j in 1:9) {
        index <- 455 * (a - 1) + 72 * 6 + 8 + j 
        if (epFlag) {
          sampleWeeklyIncidence[weekNo + 1, 9 * (a - 1) + j] <- x0[index] * ep_t[a]
        } else {
          sampleWeeklyIncidence[weekNo + 1, 9 * (a - 1) + j] <- x0[index] # Incidence at t_d = 7
        }
      }
    }
    weekNo <- weekNo + 1
  }
  
  dayNoAfterBurn <- dayNoAfterBurn + 1
  
  return(list(x0 = x0, sampleWeeklyIncidence = sampleWeeklyIncidence, dayNoAfterBurn = dayNoAfterBurn, weekNo = weekNo))
}
###
getMonthlyLikelihood <- function(x0, inc_tot, dayNoAfterBurn, monthNo, A, loglikelihoodError, valueLogLikelihood) {
  if (dayNoAfterBurn == 0) {
    x0[seq(23, by = 23, length.out = A) + 22] <- 0 # Reset incidence at t_d = 0
  }
  
  if (dayNoAfterBurn %% 30 == 0 && dayNoAfterBurn > 0) {
    for (a in 1:A) {
      index <- 22 + a * 23
      modelIncidencePerTime <- x0[index] # Extract model incidence for the month
      inc_tot[a] <- inc_tot[a] + modelIncidencePerTime
      x0[index] <- 0 # Reset for next period
    }
    # evaluateTimeStepLogLikelihood returns a list with logLikelihood and logLikelihoodError
    evaluateResult <- evaluateTimeStepLogLikelihood(monthNo) # Placeholder for actual function call
    valueLogLikelihood <- valueLogLikelihood + evaluateResult$valueLogLikelihood
    loglikelihoodError <- evaluateResult$loglikelihoodError || loglikelihoodError
    
    if (loglikelihoodError) {
      return(list(x0 = x0, inc_tot = inc_tot, dayNoAfterBurn = dayNoAfterBurn, monthNo = monthNo, loglikelihoodError = loglikelihoodError, valueLogLikelihood = valueLogLikelihood))
    }
    monthNo <- monthNo + 1
  }
  
  # Annual reset of inc_tot
  if (dayNoAfterBurn %% 360 == 0 && dayNoAfterBurn > 0) {
    # ll += check_incidence(inc_tot)
    inc_tot <- numeric(A) # Reset inc_tot for the new year
  }
  
  dayNoAfterBurn <- dayNoAfterBurn + 1
  
  return(list(x0 = x0, inc_tot = inc_tot, dayNoAfterBurn = dayNoAfterBurn, monthNo = monthNo, loglikelihoodError = loglikelihoodError, valueLogLikelihood = valueLogLikelihood))
}
###
getMonthlyIncidence <- function(x0, sampleMonthlyIncidence, epFlag, dayNoAfterBurn, monthNo, A, ep_t) {
  if (dayNoAfterBurn == 0) {
    x0[seq(23, by = 23, length.out = A) + 22] <- 0 # Reset incidence at t_d = 0
  }
  
  if (dayNoAfterBurn %% 30 == 0 && dayNoAfterBurn > 0) {
    for (a in 1:A) {
      index <- 22 + a * 23
      # Update sampleMonthlyIncidence based on epFlag
      if (epFlag) {
        sampleMonthlyIncidence[monthNo + 1, a] <- x0[index] * ep_t[a] 
      } else {
        sampleMonthlyIncidence[monthNo + 1, a] <- x0[index] # No adjustment
      }
      x0[index] <- 0 # Reset for next period
    }
    monthNo <- monthNo + 1
  }
  
  dayNoAfterBurn <- dayNoAfterBurn + 1
  
  # Return updated variables
  return(list(x0 = x0, sampleMonthlyIncidence = sampleMonthlyIncidence, dayNoAfterBurn = dayNoAfterBurn, monthNo = monthNo))
}
### 
getAnnualIncidence <- function(x0, sampleAnnualIncidence, dayNoAfterBurn, weekNo, A) {
  if (dayNoAfterBurn == 0) {
    indicesToReset <- seq(23, by = 23, length.out = A) + 22
    x0[indicesToReset] <- 0 # Reset incidence at t_d = 0
  }
  
  if (dayNoAfterBurn %% 7 == 0 && dayNoAfterBurn > 0) {
    for (a in 1:A) {
      index <- 22 + a * 23
      sampleAnnualIncidence[a] <- sampleAnnualIncidence[a] + x0[index] # Update annual incidence
      x0[index] <- 0 # Reset for next period
    }
    weekNo <- weekNo + 1
  }
  
  dayNoAfterBurn <- dayNoAfterBurn + 1
  
  return(list(x0 = x0, sampleAnnualIncidence = sampleAnnualIncidence, dayNoAfterBurn = dayNoAfterBurn, weekNo = weekNo))
}
###
create_ODE_desc <- function(finELL) {
  xi <- 1.0 / finELL$parameterValues["xi"]
  si <- 1.0 / finELL$parameterValues["si"]
  ga0 <- 1.0 / finELL$parameterValues["ga0"]
  ga1 <- 1.0 / (finELL$parameterValues["ga0"] * finELL$parameterValues["g1"])
  ga2 <- 1.0 / (finELL$parameterValues["ga0"] * finELL$parameterValues["g1"] * finELL$parameterValues["g2"])
  ga3 <- ga2
  om <- 1.0 / finELL$parameterValues["om"]
  
  # Initialize other parameters directly from finELL
  rho <- 1.0
  alpha_i <- finELL$parameterValues["alpha_i"]
  d1 <- finELL$parameterValues["d1"]
  d2 <- d1 * finELL$parameterValues["d2"]
  d3 <- d1 * finELL$parameterValues["d2"] * finELL$parameterValues["d3"]
  a1 <- a2 <- a3 <- 1.0
  
  phi <- finELL$parameterValues["phi"]
  qp <- finELL$parameterValues["qp"]
  qc <- finELL$parameterValues["qc"]
  b1 <- finELL$parameterValues["b1"]
  psi <- finELL$parameterValues["psi"]
  
  A <- finELL$A
  eta <- finELL$eta
  populationPerAgeGroup <- finELL$populationPerAgeGroup
  contactMatrixPhy <- finELL$contactMatrixPhy
  contactMatrixCon <- finELL$contactMatrixCon
  dailyBirthRate <- finELL$dailyBirthRate
  pA <- finELL$pA
  
  # Return a list that represents the ODE_desc object with all initialized variables
  return(list(xi = xi, si = si, ga0 = ga0, ga1 = ga1, ga2 = ga2, ga3 = ga3, om = om,
              rho = rho, alpha_i = alpha_i, d1 = d1, d2 = d2, d3 = d3, a1 = a1, a2 = a2, a3 = a3,
              phi = phi, qp = qp, qc = qc, b1 = b1, psi = psi, A = A, eta = eta,
              populationPerAgeGroup = populationPerAgeGroup, contactMatrixPhy = contactMatrixPhy,
              contactMatrixCon = contactMatrixCon, dailyBirthRate = dailyBirthRate, pA = pA))
}
###