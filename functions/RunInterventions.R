### source functions 
# source("/Users/vincentarumadri/Desktop/Epi/Modelling/master_thesis_ua/functions/ODE_desc.R")

####
ParameterValuesforODE <- function(currentParamValues) {
  if(length(currentParamValues) != 25) {
    stop("currentParamValues must have 25 elements.")
  }
  
  parameterValuesTemp <- numeric(25)
  
  names(parameterValuesTemp) <- c("xi", "si", "ga0", "g1", "g2", "om", "pA1", "pA2", "pA3", "pA4", 
                                  "alpha_i", "d1", "d2", "d3", "phi", "qp", "qc", "b1", "psi", 
                                  "c5ep1", "c5ep2", "ep5", "ep6", "I1", "I2")
  
  parameterValuesTemp <- setNames(currentParamValues, names(parameterValuesTemp))
  
  # Define and compute ep_t and pA
  ep_t <- numeric(25)
  pA <- numeric(25)
  parameterValues <- parameterValuesTemp
  
  for (a in 1:16) {
    ep_t[a] <- exp(parameterValues["c5ep1"] + (a - 1) * parameterValues["c5ep2"])
  }
  for (a in 17:23) {
    ep_t[a] <- parameterValues["ep5"]
  }
  for (a in 24:25) {
    ep_t[a] <- parameterValues["ep6"]
  }
  
  pA[1:12] <- rep(parameterValues["pA1"], 12)
  pA[13:16] <- rep(parameterValues["pA2"], 4)
  pA[17:18] <- rep(parameterValues["pA3"], 2)
  pA[19:25] <- rep(parameterValues["pA4"], 7)
  
  assign("ep_t", ep_t, envir = .GlobalEnv)
  assign("pA", pA, envir = .GlobalEnv)
  assign("parameterValues", parameterValues, envir = .GlobalEnv)
  
}
initial_M <- function(parameterValues, ageStratification, populationPerAgeGroup) {
  xi <- 1.0 / parameterValues[['xi']]
  
  init_con <- numeric(0)
  
  for (i in 1:(length(ageStratification) - 1)) {
    
    cdf_lower <- pexp((365 * ageStratification[i]), rate = xi)
    cdf_upper <- pexp((365 * ageStratification[i + 1]), rate = xi)
    
    init_con_temp <- (cdf_upper - cdf_lower) / ((365 * ageStratification[i + 1] - 365 * ageStratification[i]) * xi)
    init_con[i] <- init_con_temp * populationPerAgeGroup[i]
  }
  
  cdf_last <- pexp(365 * 90, rate = xi) - pexp(365 * ageStratification[length(ageStratification)], rate = xi)
  last_age_group_init_con <- cdf_last / ((365 * 90 - 365 * ageStratification[length(ageStratification)]) * xi) * populationPerAgeGroup[length(ageStratification)]
  
  init_con <- c(init_con, last_age_group_init_con)
  
  return(init_con)
}
####
poisson_cdf <- function(l, a, x) {
  if (l == 0.0 || a == 0.0) {
    return(ppois(x, lambda = 0.000001))
  } else {
    return(ppois(x, lambda = l * a))
  }
}
####
initialProportionExposure <- function(l, a1, a2) {
  
  prop <- vector("numeric", 4)
  prop[1] <- abs(poisson_cdf(l, a2, 0) - poisson_cdf(l, a1, 0)) / ((a2 - a1) * l)
  prop[2] <- abs(poisson_cdf(l, a2, 1) - poisson_cdf(l, a1, 1)) / ((a2 - a1) * l)
  prop[3] <- abs(poisson_cdf(l, a2, 2) - poisson_cdf(l, a1, 2)) / ((a2 - a1) * l)
  prop[4] <- 1 - sum(prop[1:3])
  
  return(prop)
}
####
generateInitialStates <- function(parameterValues, ageStratification, populationPerAgeGroup) {
  
  populationMatPro <- initial_M(parameterValues, ageStratification, populationPerAgeGroup)
  
  I1 <- parameterValues[["l1"]]
  I2 <- parameterValues[["l2"]]
  I3 <- 0.5
  si <- 1.0 / parameterValues[["si"]]
  g0 <- 1.0 / parameterValues[["g0"]]
  g1 <- 1.0 / (parameterValues[["g0"]] * parameterValues[["g1"]])
  g2 <- 1.0 / (parameterValues[["g0"]] * parameterValues[["g1"]] * parameterValues[["g2"]])
  d1 <- parameterValues[["d1"]]
  d2 <- parameterValues[["d1"]] * parameterValues[["d2"]]
  d3 <- parameterValues[["d1"]] * parameterValues[["d2"]] * parameterValues[["d3"]]
  
  initialStates <- c()
  for (a in 1:25) {
    if (a < 25) {
      a1 <- ageStratification[a]
      a2 <- ageStratification[a + 1]
    } else {
      a1 <- ageStratification[a]
      a2 <- 90
    }
    propEachExposureGroup <- initialProportionExposure(I3, a1, a2)  
    pI1 <- propEachExposureGroup[1]
    pI2 <- propEachExposureGroup[2]
    pI3 <- propEachExposureGroup[3]
    pI4 <- propEachExposureGroup[4]
    age_size <- populationPerAgeGroup[a] - populationMatPro[a]
    
    initialStates_i <- c(populationMatPro[a],  # M group
                         
                         pI1 * age_size * (1 - I1) * (1 - I2),  # S0
                         pI1 * age_size * I1 * si / (si + g0),  # E0
                         pI1 * age_size * I1 * g0 / (si + g0) * pA[a],  # A0
                         pI1 * age_size * I1 * g0 / (si + g0) * (1 - pA[a]),  # I0
                         pI1 * age_size * (1 - I1) * I2,  # R0,
                         0, # V0
                         
                         pI2 * age_size * (1 - d1 * I1) * (1 - I2),  # S1
                         pI2 * age_size * d1 * I1 * si / (si + g1),  # E1
                         pI2 * age_size * d1 * I1 * g1 / (si + g1) * pA[a],  # A1
                         pI2 * age_size * d1 * I1 * g1 / (si + g1) * (1 - pA[a]),  # I1
                         pI2 * age_size * (1 - d1 * I1) * I2,  # R1
                         0, # V1
                         
                         pI3*age_size*(1.0 - d2*I1)*(1.0-I2),      # S2
                         pI3*age_size*d2*I1*si/(si+g2),            # E2
                         pI3*age_size*d2*I1*g2/(si+g2)*pA[a],      # A2
                         pI3*age_size*d2*I1*g2/(si+g2)*(1-pA[a]),  # I2
                         pI3*age_size*(1.0 - d2*I1)*I2,      # R2
                         0, # V2
                         
                         pI4*age_size*(1.0 - d3*I1)*(1.0-I2),        # S3
                         pI4*age_size*d3*I1*si/(si+g2),       # E3
                         pI4*age_size*d3*I1*g2/(si+g2)*pA[a], # A3
                         pI4*age_size*d3*I1*g2/(si+g2)*(1-pA[a]),   # I3
                         pI4*age_size*(1.0 - d3*I1)*I2, # R3
                         0,    # V3
                         
                         0  # Z
    )
    initialStates <- c(initialStates, initialStates_i) # these are vaccine states which are intiially 0
  }
  return(initialStates)
}

####
getWeeklyIncidence <- function(x0, sampleWeeklyIncidence, no_doses, epFlag) {
  if (dayNoAfterBurn %% 365 == 0) {
    dayNoAfterBurn <- 0
    for (a in 1:A) {
      indices <- ag*a + sg*6 + 1:4
      x0[indices] <- 0
      x0[(ag*a + sg*6 + 8):(ag*a + sg*6 + 16)] <- 0
    }
  }
  
  if (dayNoAfterBurn %% 7 == 0 && dayNoAfterBurn > 0) {
    for (a in 1:A) {
      for (j in 1:9) {
        index <- ag*a + sg*6 + 7 + j
        if (epFlag) {
          sampleWeeklyIncidence[weekNo, 9*(a-1) + j] <- x0[index] * ep_t[a]
        } else {
          sampleWeeklyIncidence[weekNo, 9*(a-1) + j] <- x0[index]
        }
        x0[index] <- 0
      }
    }
    
    for (a in 1:A) {
      indices <- ag*a + sg*6 + 1:4
      no_doses[weekNo, 1:4] <- no_doses[weekNo, 1:4] + x0[indices]
      x0[indices] <- 0
    }
    weekNo <- weekNo + 1
  }
  
  dayNoAfterBurn <- dayNoAfterBurn + 1
  
  return(list(x0 = x0, 
              sampleWeeklyIncidence = sampleWeeklyIncidence, 
              no_doses = no_doses, 
              dayNoAfterBurn = dayNoAfterBurn, 
              weekNo = weekNo)) 
}
################ Sample ###################
Sample <- function(vac_calendar, vac_dose, cov_c, vac_info, posteriors, run_start, run_full, run_burn, dt, A, ag, sg, ep_t) {
  require(deSolve)
  
  currentODETime <- run_start
  weekNo <- 0
  dayNoAfterBurn <- 0
  
  currentParamValues <- posteriors
  ParameterValuesforODE(currentParamValues) 
  x0 <- generateInitialStates(cov_c)
  
  ODE_desc_inst <- ODE_desc(vac_calendar, vac_dose, vac_info, cov_c)
  
  sampleWeeklyIncidence <- matrix(nrow = 521, ncol = A * 9)
  no_doses <- matrix(nrow = 521, ncol = 4)
  
  while (currentODETime < (run_full + run_burn)) {
    
    ode_result <- ode(y = x0, times = c(currentODETime, currentODETime + dt), func = model_rsv, parms = ODE_desc_inst)
    x0 <- ode_result[nrow(ode_result), -1]  
    
    if (currentODETime > run_burn) {
      results <- getWeeklyIncidence(x0, sampleWeeklyIncidence, no_doses, FALSE, dayNoAfterBurn, weekNo, A, ag, sg, ep_t)
      x0 <- results$x0
      sampleWeeklyIncidence <- results$sampleWeeklyIncidence
      no_doses <- results$no_doses
      dayNoAfterBurn <- results$dayNoAfterBurn
      weekNo <- results$weekNo
    }
    
    currentODETime <- currentODETime + dt
  }
  
  return(list(inci = sampleWeeklyIncidence, doses = no_doses))
}

### 
StatesValues <- function(vac_calendar, vac_dose, cov_c, vac_info, posteriors) {
  require(deSolve)
  
  currentODETime <- run_start
  inc_tot <- numeric(A) 
  weekNo <- 0
  dayNoAfterBurn <- 0
  
  currentParamValues <- posteriors 
  ParameterValuesforODE(currentParamValues)
  
  collect_protect <- matrix(0, nrow = 365 * 12 + 365, ncol = 30 * 6)
  x0 <- generateInitialStates(cov_c) 
  x0_N <- length(x0)
  xmat <- matrix(0, nrow = 365 * 12 + 365, ncol = x0_N)
  
  ODE_desc_inst <- ODE_desc(vac_calendar, vac_dose, vac_info, cov_c)
  
  sampleWeeklyIncidence <- matrix(0, nrow = 521, ncol = A*9)
  no_doses <- matrix(0, nrow = 521, ncol = 4)
  
  while (currentODETime < (run_full + run_burn)) {
    xmat[currentODETime + 1, ] <- x0 

    run_ode <- ode(y = x0, times = c(currentODETime, currentODETime + dt), func = model_rsv, parms = ODE_desc_inst)
    x0 <- tail(run_ode[, -1], 1) # Update state based on ODE solution
    
    if (currentODETime > run_burn) {
      
      getWeeklyIncidence(x0, sampleWeeklyIncidence, no_doses, FALSE)
    }
    
    currentODETime <- currentODETime + dt
  }
  
  return(xmat)
}

############### CwX Contacts for persons who are cocooned ###############
get_cwn <- function(prop_c, s) {
  cwn_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:12) {
      for (j in 1:A) {
        cwn_e[i, j] <- nwn_p(i, j)
      }
    }
    
    for (i in 1:12) {
      for (j in 1:12) {
        cwn_e[i, j] <- nwn_p(i, j) * (1 - prop_c)
      }
    }
    
    for (i in 18:21) {
      for (j in 1:A) {
        cwn_e[i, j] <- pwn_p(i, j)
      }
    }
    
    for (i in 18:21) {
      for (j in 1:12) {
        cwn_e[i, j] <- (cnt_matrix_p[i, j] - cnt_matrix_p_h[i, j]) * (1 - prop_c)
      }
    }
    
  } else if (s == 'c') {
    for (i in 1:12) {
      for (j in 1:A) {
        cwn_e[i, j] <- nwn_c(i, j)
      }
    }
    
    for (i in 1:12) {
      for (j in 1:12) {
        cwn_e[i, j] <- nwn_c(i, j) * (1 - prop_c)
      }
    }
    
    for (i in 18:21) {
      for (j in 1:A) {
        cwn_e[i, j] <- pwn_c(i, j)
      }
    }
    
    for (i in 18:21) {
      for (j in 1:12) {
        cwn_e[i, j] <- (cnt_matrix_c[i, j] - cnt_matrix_c_h[i, j]) * (1 - prop_c)
      }
    }
    
  } else {
    stop("Error cont")
  }
  
  return(cwn_e)
}
####
get_cwp <- function(prop_c, s) {
  cwp_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:12) {
      for (j in 18:21) {
        cwp_e[i, j] <- (cnt_matrix_p[i, j] - cnt_matrix_p_h[i, j]) * p_mat[j] * (1 - prop_c)
      }
    }
    
    for (i in 18:21) {
      for (j in 18:21) {
        cwp_e[i, j] <- pwp_p[i, j] * (1 - prop_c)
      }
    }
    
  } else if (s == 'c') {
    for (i in 1:12) {
      for (j in 18:21) {
        cwp_e[i, j] <- (cnt_matrix_c[i, j] - cnt_matrix_c_h[i, j]) * p_mat[j] * (1 - prop_c)
      }
    }
    
    for (i in 18:21) {
      for (j in 18:21) {
        cwp_e[i, j] <- pwp_c[i, j] * (1 - prop_c)
      }
    }
    
  } else {
    stop("Error cont")
  }
  
  return(cwp_e)
}
####
get_cwc <- function(prop_c, s) {
  cwc_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:12) {
      for (j in 1:12) {
        cwc_e[i, j] <- nwn_p(i, j) * prop_c
      }
    }
    
    for (i in 18:21) {
      for (j in 1:12) {
        cwc_e[i, j] <- cnt_matrix_p_h(i, j) + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j)) * prop_c
      }
    }
    
    for (i in 1:12) {
      for (j in 18:21) {
        cwc_e[i, j] <- cnt_matrix_p_h(i, j) / 2.0 + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j)) * p_mat[j] * prop_c
      }
    }
    
    for (i in 18:21) {
      for (j in 18:21) {
        cwc_e[i, j] <- pwp_p(i, j) * prop_c
      }
    }
    
  } else if (s == 'c') {
    for (i in 1:12) {
      for (j in 1:12) {
        cwc_e[i, j] <- nwn_c(i, j) * prop_c
      }
    }
    
    for (i in 18:21) {
      for (j in 1:12) {
        cwc_e[i, j] <- cnt_matrix_c_h(i, j) + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j)) * prop_c
      }
    }
    
    for (i in 1:12) {
      for (j in 18:21) {
        cwc_e[i, j] <- cnt_matrix_c_h(i, j) / 2.0 + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j)) * p_mat[j] * prop_c
      }
    }
    
    for (i in 18:21) {
      for (j in 18:21) {
        cwc_e[i, j] <- pwp_c(i, j) * prop_c
      }
    }
    
  } else {
    stop("Error cont")
  }
  
  return(cwc_e)
}
###### PwX Contacts for persons who are mothers but not cocooned  #######
get_pwn <- function(prop_c, s) {
  pwn_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 18:21) {
      for (j in 1:A) {
        pwn_e[i, j] <- pwn_p(i, j)
      }
    }
    
    for (i in 18:21) {
      for (j in 1:12) {
        pwn_e[i, j] <- cnt_matrix_p_h(i, j) + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j)) * (1 - prop_c)
      }
    }
    
  } else if (s == 'c') {
    for (i in 18:21) {
      for (j in 1:A) {
        pwn_e[i, j] <- pwn_c(i, j)
      }
    }
    
    for (i in 18:21) {
      for (j in 1:12) {
        pwn_e[i, j] <- cnt_matrix_c_h(i, j) + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j)) * (1 - prop_c)
      }
    }
    
  } else {
    stop("Error cont")
  }
  
  return(pwn_e)
}
####
get_pwp <- function(prop_c, s) {
  pwp_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 18:21) {
      for (j in 18:21) {
        pwp_e[i, j] <- pwp_p(i, j) * (1 - prop_c)
      }
    }
    
  } else if (s == 'c') {
    for (i in 18:21) {
      for (j in 18:21) {
        pwp_e[i, j] <- pwp_c(i, j) * (1 - prop_c)
      }
    }
    
  } else {
    stop("Error cont")
  }
  
  return(pwp_e)
}
####
get_pwc <- function(prop_c, s) {
  pwc_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 18:21) {
      for (j in 1:12) {
        pwc_e[i, j] <- (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j)) * prop_c
      }
    }
    
    for (i in 18:21) {
      for (j in 18:21) {
        pwc_e[i, j] <- pwn_p(i, j) * prop_c
      }
    }
    
  } else if (s == 'c') {
    for (i in 18:21) {
      for (j in 1:12) {
        pwc_e[i, j] <- (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j)) * prop_c
      }
    }
    
    for (i in 18:21) {
      for (j in 18:21) {
        pwc_e[i, j] <- pwn_c(i, j) * prop_c
      }
    }
    
  } else {
    stop("Error cont")
  }
  
  return(pwc_e)
}
##### NwX Contacts for persons who are neither cocooned nor mothers #####
get_nwn <- function(prop_c, s) {
  nwn_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:A) {
      for (j in 1:A) {
        nwn_e[i, j] <- nwn_p(i, j)
      }
      for (j in 1:12) {
        nwn_e[i, j] <- nwn_p(i, j) * (1 - prop_c)
      }
    }
  } else if (s == 'c') {
    for (i in 1:A) {
      for (j in 1:A) {
        nwn_e[i, j] <- nwn_c(i, j)
      }
      for (j in 1:12) {
        nwn_e[i, j] <- nwn_c(i, j) * (1 - prop_c)
      }
    }
  } else {
    stop("Error cont")
  }
  
  return(nwn_e)
}
####
get_nwp <- function(prop_c, s) {
  nwp_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:12) {
      for (j in 18:21) {
        nwp_e[i, j] <- (cnt_matrix_p_h[i, j] / 2.0) + (cnt_matrix_p[i, j] - cnt_matrix_p_h[i, j]) * p_mat[j] * (1 - prop_c)
      }
    }
    for (i in 13:25) {
      for (j in 18:21) {
        nwp_e[i, j] <- cnt_matrix_p[i, j] * p_mat[j] * (1 - prop_c)
      }
    }
  } else if (s == 'c') {
    for (i in 1:12) {
      for (j in 18:21) {
        nwp_e[i, j] <- (cnt_matrix_c_h[i, j] / 2.0) + (cnt_matrix_c[i, j] - cnt_matrix_c_h[i, j]) * p_mat[j] * (1 - prop_c)
      }
    }
    for (i in 13:25) {
      for (j in 18:21) {
        nwp_e[i, j] <- cnt_matrix_c[i, j] * p_mat[j] * (1 - prop_c)
      }
    }
  } else {
    stop("Error cont")
  }
  
  return(nwp_e)
}
####
get_nwc <- function(prop_c, s) {
  nwc_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:A) {
      for (j in 1:12) {
        nwc_e[i, j] <- nwn_p(i, j) * prop_c
      }
    }
    for (i in 1:12) {
      for (j in 18:21) {
        nwc_e[i, j] <- (cnt_matrix_p[i, j] - cnt_matrix_p_h[i, j]) * p_mat[j] * prop_c
      }
    }
    for (i in 13:25) {
      for (j in 18:21) {
        nwc_e[i, j] <- cnt_matrix_p[i, j] * p_mat[j] * prop_c
      }
    }
  } else if (s == 'c') {
    for (i in 1:A) {
      for (j in 1:12) {
        nwc_e[i, j] <- nwn_c(i, j) * prop_c
      }
    }
    for (i in 1:12) {
      for (j in 18:21) {
        nwc_e[i, j] <- (cnt_matrix_c[i, j] - cnt_matrix_c_h[i, j]) * p_mat[j] * prop_c
      }
    }
    for (i in 13:25) {
      for (j in 18:21) {
        nwc_e[i, j] <- cnt_matrix_c[i, j] * p_mat[j] * prop_c
      }
    }
  } else {
    stop("Error cont")
  }
  
  return(nwc_e)
}

populationPerAgeGroup <- uk_data_sum$populationAgeGroup
eta <- numeric()
################ RunInterventions ###################
RunInterventions_func <- function(dailyBirthRate, totPopulation, ageStratification) {
  dailyBirthRate <- dailyBirthRate
  totPopulation <- totPopulation
  ageStratification <- ageStratification
  
  A <- length(ageStratification)
  eta <- c(0) 
  populationPerAgeGroup <- c() 
  modelIncidencePerTime <- c() 
  
  for (i in 1:(A-1)) {
    populationPerAgeGroup <- c(populationPerAgeGroup, dailyBirthRate*365*(ageStratification[i+1] - ageStratification[i]))
    eta <- c(eta, 1.0/(365.0*(ageStratification[i+1] - ageStratification[i])))
    modelIncidencePerTime <- c(modelIncidencePerTime, 0)
  }
  
  modelIncidencePerTime <- c(modelIncidencePerTime) 
  populationPerAgeGroup <- c(populationPerAgeGroup, totPopulation - (dailyBirthRate*365)*ageStratification[A])
  eta <- c(eta, dailyBirthRate / (totPopulation - (dailyBirthRate*365)*ageStratification[A]))
  
  dt <- 1
  currentODETime <- 0
  dayNoAfterBurn <- 0
  valueLogLikelihood <- 0
  
  rg <- 45
  sg <- rg * 3
  ag <- sg * 6 + 23
  
  return(list(
    dailyBirthRate = dailyBirthRate,
    totPopulation = totPopulation,
    ageStratification = ageStratification,
    A = A,
    eta = eta,
    populationPerAgeGroup = populationPerAgeGroup,
    modelIncidencePerTime = modelIncidencePerTime,
    dt = dt,
    currentODETime = currentODETime,
    dayNoAfterBurn = dayNoAfterBurn,
    valueLogLikelihood = valueLogLikelihood,
    rg = rg,
    sg = sg,
    ag = ag
  ))
}

### foi #####
calculate_foi <- function(t, state_initial_age, contactMatrixPhy, contactMatrixCon, qp, b1, phi, psi, delta_i, qc) {
  # initialize for the current time point
  lambda_matrix <- matrix(0, nrow = 5, ncol = 4)
  contact_rate <- matrix(0, nrow = 5, ncol = 4)
  
  # initialize array to store N_b values
  N_b <- numeric(5)
  
  # each infection level 'i' and age group 'a'
  for (i in 1:4) {
    for (a in 1:5) {
      # initialize the contact matrix the current age group and infection level
      summation <- 0
      
      # N_b and the contact matrix for each 'b'
      for (b in 1:5) {
        # Calculate N_b for current 'b'
        N_bi <- 0
        for (c in 0:3) {
          S_name <- paste0("S", c, b)
          matching_elements <- state_initial_age[grepl(S_name, names(state_initial_age))]
          if (length(matching_elements) > 0) {
            N_bi <- N_bi + sum(matching_elements)
          }
        }
        N_b[b] <- N_bi  # store the computed N_b
        
        # retrieve the A and I values for current infection level 'i' and age group 'b'
        A_bi <- state_initial_age[paste0("A", i-1, b)]
        I_bi <- state_initial_age[paste0("I", i-1, b)]
        
        # Add the contribution of age group 'b' to the summation
        if (N_b[b] > 0) {  # avoid division by zero
          contacts <- ((contactMatrixPhy[a, b]) + (contactMatrixCon[a, b] * qc)) * (A_bi + I_bi) / N_b[b]
          summation <- summation + contacts
        }
      }
      
      contact_rate[a, i] <- summation
      
      # product of delta terms up to the infection level 'i'
      delta_product <- prod(delta_i[1:i])
      
      t1 = t %% 365
      # Calculate the time-varying factor
      beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))
      
      # Calculate lambda for the current age group 'a' and infection level 'i'
      lambda_matrix[a, i] <- qp * beta_value * delta_product * summation
    }
  }
  
  return(lambda_matrix)
}

calculate_gamma <- function(params) {
  
  # parameters from list
  g0 = params[['g0']]
  g1 = params[['g1']]
  g2 = params[['g2']]
  
  # new parameters
  params[['gamma0']] <<- 1/g0
  params[['gamma1']] <<- 1/(g0*g1)
  params[['gamma2']] <<- 1/(g0*g1*g2)
  
}
plot_figures <- function(ode_output, programme_name){
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse")
    library(tidyverse)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
    library(patchwork)
  }
  
  ode_output <- as.data.frame(ode_output)
  programme_name <- as.character(programme_name)
  
  # p_maternal
  p_maternal <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("M1","M2", "M3"))  
  
  p_maternal <- ggplot(p_maternal, aes(x = time, y = value, color = state)) +
    geom_line() +
    labs(title = "Maternal",
         x = "Time",
         y = "Population") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,1825), breaks = c(365,730,1095,1460, 1825),
                       labels = c("1yr", "2", "3", "4", "5") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"), 
          plot.tag=element_text(face="bold"),
          legend.position = "none") 
  # level 0
  p_level0 <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("S01", "S02", "S03", "E01", "E02", "E03", "A01", "A02", "A03", 
                        "I01", "I02", "I03", "R01", "R02", "R03", "V01", "V02", "V03")) %>%
    mutate(state = factor(state, levels = c("S01", "S02", "S03", "E01", "E02", "E03", "A01", "A02", "A03",
                                            "I01", "I02", "I03", "R01", "R02", "R03", "V01", "V02", "V03")))
  
  p_level0 <- ggplot(p_level0, aes(x = time, y = value, color = state)) +
    geom_line() +
    labs(title = "Exposure level 0",
         x ="",
         y = "Population",
         tag = "A") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,1825), breaks = c(365,730,1095,1460, 1825),
                       labels = c("1yr", "2", "3", "4", "5") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"),
          plot.tag=element_text(face="bold"),
          legend.position = "none") 
  
  # level 1
  p_level1 <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("S11", "S12", "S13", "E11", "E12", "E13", "A11", "A12", "A13", 
                        "I11", "I12", "I13", "R11", "R12", "R13", "V11", "V12", "V13"))%>%
    mutate(state = factor(state, levels = c("S11", "S12", "S13", "E11", "E12", "E13", "A11", "A12", "A13",
                                            "I11", "I12", "I13", "R11", "R12", "R13", "V11", "V12", "V13")))
  
  p_level1 <- ggplot(p_level1, aes(x = time, y = value, color = state)) +
    geom_line() +
    labs(title = "Exposure level 1",
         x = "",
         y = "",
         tag = "B") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,1825), breaks = c(365,730,1095,1460, 1825),
                       labels = c("1yr", "2", "3", "4", "5") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"),
          plot.tag=element_text(face="bold"),
          legend.position = "none") 
  
  # level 2
  p_level2 <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("S21", "S22", "S23", "E21", "E22", "E23", "A21", "A22", "A23", 
                        "I21", "I22", "I23", "R21", "R22", "R23", "V21", "V22", "V23"))%>%
    mutate(state = factor(state, levels = c("S21", "S22", "S23", "E21", "E22", "E23", "A21", "A22", "A23",
                                            "I21", "I22", "I23", "R21", "R22", "R23", "V21", "V22", "V23")))
  
  p_level2 <- ggplot(p_level2, aes(x = time, y = value, color = state)) +
    geom_line() +
    labs(title = "Exposure level 2",
         x = "time",
         y = "Population",
         tag = "C") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,1825), breaks = c(365,730,1095,1460, 1825),
                       labels = c("1yr", "2", "3", "4", "5") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"),
          plot.tag=element_text(face="bold"),
          legend.position = "none") 
  
  # level 3
  p_level3 <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("S31", "S32", "S33", "E31", "E32", "E33", "A31", "A32", "A33", 
                        "I31", "I32", "I33", "R31", "R32", "R33", "V31", "V32", "V33"))%>%
    mutate(state = factor(state, levels = c("S31", "S32", "S33", "E31", "E32", "E33", "A31", "A32", "A33",
                                            "I31", "I32", "I33", "R31", "R32", "R33", "V31", "V32", "V33")))
  
  p_level3 <- ggplot(p_level3, aes(x = time, y = value, color = state)) +
    geom_line() +
    labs(title = "Exposure level 3",
         x = "time",
         y = "",
         tag = "D") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,1825), breaks = c(365,730,1095,1460, 1825),
                       labels = c("1yr", "2", "3", "4", "5") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"),
          plot.tag=element_text(face="bold"),
          legend.position = "none" ) 
  
  # joint 
  plot_base <- (p_level0 | p_level1) / (p_level2 | p_level3)
  
  plot_base <- plot_base + plot_annotation(
    title = programme_name,
    subtitle = "Five year movement of individuals between disease states and across age groups",
    theme = theme(
      plot.title = element_text(family="sans",size=15, color = "black"),
      plot.subtitle = element_text(family="sans",size=12, color = "black")+
        theme_classic()
    )
  )
  
  return(plot_base)
}

new_cases <- function(ode_output){
  start_simulation <- (6*365 + 1)
  end_simulation <- 3650
  
  relevant_data <- ode_output[start_simulation:end_simulation, ]
  
  daily_changes <- diff(relevant_data[,c("Z1", "Z2", "Z3")])
  daily_changes_df <- as.data.frame(daily_changes)
  
  new_cases <- colSums(daily_changes)
  
  return(new_cases)
}

annual_incidence <- function(ode_output){
  start_simulation <- (8*365 + 1)
  end_simulation <- 3285
  
  relevant_data <- ode_output[start_simulation:end_simulation, ]
  
  daily_changes <- diff(relevant_data[,c("Z1", "Z2", "Z3")])
  daily_changes_df <- as.data.frame(daily_changes)
  
  new_cases <- colSums(daily_changes)
  
  return(new_cases)
}
