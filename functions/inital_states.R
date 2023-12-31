####
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
  }
  return(initialStates)
}
