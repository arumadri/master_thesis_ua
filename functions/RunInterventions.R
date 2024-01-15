############### CwX Contacts for persons who are cocooned ###############
get_cwn <- function(prop_c, s, A, cnt_matrix_p, cnt_matrix_p_h, nwn_p, pwn_p, cnt_matrix_c, nwn_c, pwn_c) {
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
get_cwp <- function(prop_c, s, A, cnt_matrix_p, cnt_matrix_p_h, p_mat, pwp_p, cnt_matrix_c, pwp_c) {
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
get_cwc <- function(prop_c, s, A, cnt_matrix_p, cnt_matrix_p_h, nwn_p, pwp_p, cnt_matrix_c, nwn_c, pwp_c, p_mat) {
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
get_pwn <- function(prop_c, s, A, cnt_matrix_p, cnt_matrix_p_h, pwn_p, cnt_matrix_c, pwn_c) {
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
get_pwp <- function(prop_c, s, A, pwp_p, pwp_c) {
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
get_pwc <- function(prop_c, s, A, cnt_matrix_p, cnt_matrix_p_h, pwn_p, cnt_matrix_c, cnt_matrix_c_h, pwn_c) {
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
get_nwn <- function(prop_c, s, A, nwn_p, nwn_c) {
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
get_nwp <- function(prop_c, s, A, cnt_matrix_p, cnt_matrix_p_h, cnt_matrix_c, cnt_matrix_c_h, p_mat) {
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
get_nwc <- function(prop_c, s, A, nwn_p, cnt_matrix_p, cnt_matrix_p_h, p_mat, nwn_c, cnt_matrix_c, cnt_matrix_c_h) {
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
# Initialize vectors
populationPerAgeGroup <- numeric()
eta <- numeric()
modelIncidencePerTime <- numeric()
pA <- numeric()
ep_t <- numeric()

# Initialize numeric variables
dailyBirthRate <- 0
totPopulation <- 0
valueLogLikelihood <- 0

# Defined later but of size A
currentODETime <- 0
run_start <- 0
run_burn <- 0
run_full <- 0
dt <- 0

# Initialize ageStratification as a numeric vector
ageStratification <- numeric()

# Initialize integer variables
dayNoAfterBurn <- 0
weekNo <- 0
monthNo <- 0
################ RunInterventions ###################
RunInterventions <- function(dailyBirthRate_t, totPopulation_t, ageStratification_t) {
  # variables 
  dailyBirthRate <- dailyBirthRate_t
  totPopulation <- totPopulation_t
  ageStratification <- ageStratification_t
  A <- length(ageStratification)
  eta <- numeric(A)
  populationPerAgeGroup <- numeric(A)
  modelIncidencePerTime <- numeric(A)
  eta[1] <- 0
  
  # Using for loop for calculations
  for (i in 1:(A - 1)) {
    populationPerAgeGroup[i] <- dailyBirthRate * 365 * (ageStratification[i + 1] - ageStratification[i])
    eta[i + 1] <- 1.0 / (365.0 * (ageStratification[i + 1] - ageStratification[i]))
  }
  modelIncidencePerTime[1:(A - 1)] <- 0
  populationPerAgeGroup[A] <- totPopulation - (dailyBirthRate * 365) * ageStratification[A]
  eta[A] <- dailyBirthRate / (totPopulation - (dailyBirthRate * 365) * ageStratification[A])
  modelIncidencePerTime[A] <- 0
  
  # Other variables 
  dt <- 1
  currentODETime <- 0
  dayNoAfterBurn <- 0
  valueLogLikelihood <- 0
  
  # new values
  list(
    A = A,
    eta = eta,
    populationPerAgeGroup = populationPerAgeGroup,
    modelIncidencePerTime = modelIncidencePerTime,
    dt = dt,
    currentODETime = currentODETime,
    dayNoAfterBurn = dayNoAfterBurn,
    valueLogLikelihood = valueLogLikelihood
  )
}

# initialise pVHR, pHR, pLR
pVHR <- numeric() 
pHR <- numeric()   
pLR <- numeric()
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
  
  for (a in 1:16) {
    ep_t[a] <- exp(parameterValuesTemp["c5ep1"] + (a - 1) * parameterValuesTemp["c5ep2"])
  }
  for (a in 17:23) {
    ep_t[a] <- parameterValuesTemp["ep5"]
  }
  for (a in 24:25) {
    ep_t[a] <- parameterValuesTemp["ep6"]
  }
  
  pA[1:12] <- rep(parameterValuesTemp["pA1"], 12)
  pA[13:16] <- rep(parameterValuesTemp["pA2"], 4)
  pA[17:18] <- rep(parameterValuesTemp["pA3"], 2)
  pA[19:25] <- rep(parameterValuesTemp["pA4"], 7)
  
  list(ep_t = ep_t, pA = pA, parameterValues = parameterValuesTemp)
}
initial_M <- function(parameterValues, ageStratification, populationPerAgeGroup) {
  xi <- 1.0 / parameters["xi"]
  
  init_con <- numeric(length(ageStratification) - 1)
  
  for (i in 1:(length(ageStratification) - 1)) {
    cdf_upper <- pexp(365 * ageStratification[i + 1], rate = xi)
    cdf_lower <- pexp(365 * ageStratification[i], rate = xi)
    init_con_temp <- (cdf_upper - cdf_lower) / ((365 * ageStratification[i + 1] - 365 * ageStratification[i]) * xi)
    init_con[i] <- init_con_temp * populationPerAgeGroup[i]
  }
  
  cdf_upper_last <- pexp(365 * 90, rate = xi)
  cdf_lower_last <- pexp(365 * ageStratification[length(ageStratification)], rate = xi)
  init_con[length(init_con)] <- (cdf_upper_last - cdf_lower_last) / ((365 * 90 - 365 * ageStratification[length(ageStratification)]) * xi) * populationPerAgeGroup[length(populationPerAgeGroup)]
  
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
initialProportionExposure <- function(l, a1, a2, A) {
  # Assuming poisson_cdf function is already defined
  prop <- numeric(A)
  
  prop[1] <- abs(poisson_cdf(l, a2, 0) - poisson_cdf(l, a1, 0)) / ((a2 - a1) * l)
  prop[2] <- abs(poisson_cdf(l, a2, 1) - poisson_cdf(l, a1, 1)) / ((a2 - a1) * l)
  prop[3] <- abs(poisson_cdf(l, a2, 2) - poisson_cdf(l, a1, 2)) / ((a2 - a1) * l)
  prop[4] <- 1 - (prop[3] + prop[2] + prop[1])
  
  return(prop)
}
####
generateInitialStates <- function(cov_c, A, ageStratification, parameters, populationPerAgeGroup, p_mat, pVHR, pHR, pLR) {
  populationMatPro <- initial_M(parameters, ageStratification, populationPerAgeGroup)  
  initialStates <- numeric()
  
  I1 <- parameters["I1"]
  I2 <- parameters["I2"]
  I3 <- 0.5
  si <- 1.0 / parameters["si"]
  g0 <- 1.0 / parameters["ga0"]
  g1 <- 1.0 / (parameters["ga0"] * parameters["g1"])
  g2 <- 1.0 / (parameters["ga0"] * parameters["g1"] * parameters["g2"])
  d1 <- parameters["d1"]
  d2 <- parameters["d1"] * parameters["d2"]
  d3 <- parameters["d1"] * parameters["d2"] * parameters["d3"]
  
  for (a in 1:A) {
    if (a < A) {
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
    initialStates <- c(initialStates, rep(0, 23)) # inf exposure 1
  }
  return(initialStates)
}
####
getWeeklyIncidence <- function(x0, sampleWeeklyIncidence, no_doses, epFlag, A, ep_t, dayNoAfterBurn, weekNo) {
  
  if (dayNoAfterBurn == 0) {
    for (a in 1:A) {
      # ensure r indexing
      baseIndex <- 455 * (a - 1) + 72 * 6
      # Resetting specific elements of x0 to 0
      x0[baseIndex + 2] <- 0
      x0[baseIndex + 3] <- 0
      x0[baseIndex + 4] <- 0
      x0[baseIndex + 5] <- 0
      
      # Resetting a range of elements in x0 to 0
      for (j in 1:9) {
        x0[baseIndex + 8 + j] <- 0
      }
    }
  }
  
  if (dayNoAfterBurn %% 7 == 0 && dayNoAfterBurn > 0) {
    for (a in 1:A) {
      baseIndex <- 455 * (a - 1) + 72 * 6
      for (j in 1:9) {
        incidenceIndex <- baseIndex + 8 + j
        if (epFlag) {
          sampleWeeklyIncidence[weekNo, 9 * (a - 1) + j] <- x0[incidenceIndex] * ep_t[a]
        } else {
          sampleWeeklyIncidence[weekNo, 9 * (a - 1) + j] <- x0[incidenceIndex]
        }
        x0[incidenceIndex] <- 0
      }
      
      doseIndices <- (baseIndex + 2):(baseIndex + 5)
      no_doses[weekNo, 1:4] <- no_doses[weekNo, 1:4] + x0[doseIndices]
      x0[doseIndices] <- 0
    }
    
    weekNo <- weekNo + 1
  }
  
  dayNoAfterBurn <- dayNoAfterBurn + 1
  
  # Return the updated values
  list(x0 = x0, sampleWeeklyIncidence = sampleWeeklyIncidence, no_doses = no_doses, dayNoAfterBurn = dayNoAfterBurn, weekNo = weekNo)
}
################ Sample ###################
####
Sample <- function(vac_calendar, vac_dose, cov_c, vac_info, posteriors) {
  
  # Restart values
  currentODETime <- run_start
  inc_tot <- numeric(A)
  weekNo <- 0
  dayNoAfterBurn <- 0
  
  # Assign parameter values
  currentParamValues <- posteriors
  ParameterValuesforODE(currentParamValues)
  
  # Set up initial states and ODE parameters
  x0 <- generateInitialStates(cov_c, )
  ODE_desc_inst <- list(vac_calendar = vac_calendar, vac_dose = vac_dose, vac_info = vac_info, cov_c = cov_c)
  
  sampleWeeklyIncidence <- matrix(nrow = 521, ncol = A * 9)
  no_doses <- matrix(nrow = 521, ncol = 4)
  
  # Run ODE solver (conceptual, using deSolve or similar package)
  while (currentODETime < (run_full + run_burn)) {
    
    # Update x0 using an ODE solver function (e.g., ode from deSolve)
    # This is a placeholder - actual implementation depends on the ODE system
    x0 <- ode(y = x0, times = c(currentODETime, currentODETime + dt), func = ODESystemFunction, parms = ODE_desc_inst)
    
    if (currentODETime > run_burn) {
      results <- getWeeklyIncidence(x0, sampleWeeklyIncidence, no_doses, FALSE, )
      x0 <- results$x0
      sampleWeeklyIncidence <- results$sampleWeeklyIncidence
      no_doses <- results$no_doses
    }
    
    currentODETime <- currentODETime + dt
  }
  
  list(inci = sampleWeeklyIncidence, doses = no_doses)
}