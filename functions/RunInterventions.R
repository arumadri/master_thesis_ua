### source functions 
# source("/Users/vincentarumadri/Desktop/Epi/Modelling/master_thesis_ua/functions/ODE_desc.R")

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
  prop <- numeric()
  
  prop <- vector("numeric", 4)
  prop[1] <- abs(poisson_cdf(l, a2, 0) - poisson_cdf(l, a1, 0)) / ((a2 - a1) * l)
  prop[2] <- abs(poisson_cdf(l, a2, 1) - poisson_cdf(l, a1, 1)) / ((a2 - a1) * l)
  prop[3] <- abs(poisson_cdf(l, a2, 2) - poisson_cdf(l, a1, 2)) / ((a2 - a1) * l)
  prop[4] <- 1 - sum(prop[1:3])
  
  return(prop)
}
####
generateInitialStates <- function(cov_c) {
  populationMatPro <- initial_M(parameterValues, ageStratification, populationPerAgeGroup)
  initialStates <- c()
  
  I1 <- parameterValues[["I1"]]
  I2 <- parameterValues[["I2"]]
  I3 <- 0.5
  si <- 1.0 / parameterValues[["si"]]
  g0 <- 1.0 / parameterValues[["ga0"]]
  g1 <- 1.0 / (parameterValues[["ga0"]] * parameterValues[["g1"]])
  g2 <- 1.0 / (parameterValues[["ga0"]] * parameterValues[["g1"]] * parameterValues[["g2"]])
  d1 <- parameterValues[["d1"]]
  d2 <- parameterValues[["d1"]] * parameterValues[["d2"]]
  d3 <- parameterValues[["d1"]] * parameterValues[["d2"]] * parameterValues[["d3"]]
  
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
    
    for (s in 0:5) {
      if (a < 12) {
        if (s == 0 || s == 3) {
          s_prop <- 0.0
        } else if (s == 1 || s == 4) {
          s_prop <- cov_c
        } else if (s == 2 || s == 5) {
          s_prop <- 1 - cov_c
        } else {
          stop("OOPS")
        }
      } else {
        if (s == 0 || s == 3) {
          s_prop <- p_mat[a] * (1 - cov_c)
        } else if (s == 1 || s == 4) {
          s_prop <- p_mat[a] * cov_c
        } else if (s == 2 || s == 5) {
          s_prop <- 1 - p_mat[a]
        } else {
          stop("OOPS")
        }
      }
      
      for (r in 0:2) {
        if (r == 0) {
          r_prop <- pVHR[a]
        } else if (r == 1) {
          r_prop <- pHR[a]
        } else if (r == 2) {
          r_prop <- pLR[a]
        } else {
          stop("ERROR")
        }
        
        initialStates <- c(r_prop * s_prop * populationMatPro[a],  # M group
                           r_prop * s_prop * pI1 * age_size * (1 - I1) * (1 - I2),  # Sus
                           r_prop * s_prop * pI1 * age_size * I1 * si / (si + g0),  # Exp
                           r_prop * s_prop * pI1 * age_size * I1 * g0 / (si + g0) * pA[a],  # Inf A
                           r_prop * s_prop * pI1 * age_size * I1 * g0 / (si + g0) * (1 - pA[a]),  # Inf S
                           r_prop * s_prop * pI1 * age_size * (1 - I1) * I2,  # Rec
                           
                           r_prop * s_prop * pI2 * age_size * (1 - d1 * I1) * (1 - I2),  # Sus
                           r_prop * s_prop * pI2 * age_size * d1 * I1 * si / (si + g1),  # Exp
                           r_prop * s_prop * pI2 * age_size * d1 * I1 * g1 / (si + g1) * pA[a],  # Inf A
                           r_prop * s_prop * pI2 * age_size * d1 * I1 * g1 / (si + g1) * (1 - pA[a]),  # Inf S
                           r_prop * s_prop * pI2 * age_size * (1 - d1 * I1) * I2,  # Rec
                           
                           r_prop*s_prop*pI3*age_size*(1.0 - d2*I1)*(1.0-I2),      # S
                           r_prop*s_prop*pI3*age_size*d2*I1*si/(si+g2),            # Exp
                           r_prop*s_prop*pI3*age_size*d2*I1*g2/(si+g2)*pA[a],      # Inf A
                           r_prop*s_prop*pI3*age_size*d2*I1*g2/(si+g2)*(1-pA[a]),  # Inf S
                           r_prop*s_prop*pI3*age_size*(1.0 - d2*I1)*I2,      # Rec
                           
                           r_prop*s_prop*pI4*age_size*(1.0 - d3*I1)*(1.0-I2),        # Sus
                           r_prop*s_prop*pI4*age_size*d3*I1*si/(si+g2),       # Exp
                           r_prop*s_prop*pI4*age_size*d3*I1*g2/(si+g2)*pA[a], # Inf A
                           r_prop*s_prop*pI4*age_size*d3*I1*g2/(si+g2)*(1-pA[a]),   # Inf S
                           r_prop*s_prop*pI4*age_size*(1.0 - d3*I1)*I2,     # Rec
                           
                           0, # Rec
                           0, # Rec
                           0  # Rec
        )
      }
    }
    #initialStates <- c(initialStates, rep(0, 23)) # inf exposure 1
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

