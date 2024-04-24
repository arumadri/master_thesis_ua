
generateInitialStates <- function(cov_c, parameterValues, ageStratification, populationPerAgeGroup) {

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
    
    for (s in 0:5) {
      if (a < 12) {
        if (s == 0 || s == 3) {
          s_prop <- 0.0
        } else if (s == 1 || s == 4) {
          s_prop <- cov_c[a]
        } else if (s == 2 || s == 5) {
          s_prop <- 1 - cov_c[a]
        } else {
          stop("OOPS")
        }
      } else {
        if (s == 0 || s == 3) {
          s_prop <- p_mat[a] * (1 - cov_c[a])
        } else if (s == 1 || s == 4) {
          s_prop <- p_mat[a] * cov_c[a]
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
          stop("ERROR")}
        
        
        initialStates_i <- c(r_prop * s_prop * populationMatPro[a],  # M group
                           
                             r_prop * s_prop * pI1 * age_size * (1 - I1) * (1 - I2),  # S
                             r_prop * s_prop * pI1 * age_size * I1 * si / (si + g0),  # E
                             r_prop * s_prop * pI1 * age_size * I1 * g0 / (si + g0) * pA[a],  # A
                             r_prop * s_prop * pI1 * age_size * I1 * g0 / (si + g0) * (1 - pA[a]),  # I
                             r_prop * s_prop * pI1 * age_size * (1 - I1) * I2,  # R
                           
                             r_prop * s_prop * pI2 * age_size * (1 - d1 * I1) * (1 - I2),  # S
                             r_prop * s_prop * pI2 * age_size * d1 * I1 * si / (si + g1),  # E
                             r_prop * s_prop * pI2 * age_size * d1 * I1 * g1 / (si + g1) * pA[a],  # A
                             r_prop * s_prop * pI2 * age_size * d1 * I1 * g1 / (si + g1) * (1 - pA[a]),  # I
                             r_prop * s_prop * pI2 * age_size * (1 - d1 * I1) * I2,  # R
                           
                             r_prop * s_prop * pI3 * age_size * (1 - d2 * I1) * (1.0-I2),      # S
                             r_prop * s_prop * pI3 * age_size * d2 * I1 * si / (si+g2),            # E
                             r_prop * s_prop * pI3 * age_size * d2 * I1 * g2 / (si+g2) * pA[a],      # A
                             r_prop * s_prop * pI3 * age_size * d2 * I1 * g2 / (si+g2) * (1-pA[a]),  # I
                             r_prop * s_prop * pI3 * age_size * (1 - d2 * I1) * I2,      # R
                           
                             r_prop * s_prop * pI4 * age_size * (1 - d3 * I1) * (1.0-I2),        # S
                             r_prop * s_prop * pI4 * age_size * d3 * I1 * si / (si+g2),       # E
                             r_prop * s_prop * pI4 * age_size * d3 * I1 * g2 / (si+g2) * pA[a], # A
                             r_prop * s_prop * pI4 * age_size * d3 * I1 * g2 / (si+g2) * (1-pA[a]),   # I
                             r_prop * s_prop * pI4 * age_size * (1 - d3 * I1) * I2    # R
        )
         initialStates <- c(initialStates, initialStates_i, rep(0, 24)) # these are vaccine states which are intially 0
      }
    }
    initialStates <- c(initialStates, rep(0, 23)) # these are summary states for output
  }
  return(initialStates)
}



pVHR <- uk_data_sum$pVHR
pHR <- uk_data_sum$pHR
pLR <- uk_data_sum$pLR
p_mat <- uk_data_sum$prop_mat
cov_c <- uk_data_sum$nmat
initialStates_i <- generateInitialStates(p_mat, param_means, ageStratification, populationPerAgeGroup)
