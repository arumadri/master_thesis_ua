# parameters and initial condition variables 
############### parameters ##########################
parameters <- list(b1 = 1.545721,                           # mean seasonal transmission rate
                   phi = 0.6924976,                         # seasonal offset/peak transmission seasonal
                   psi = 5.068942,                          # width of seasonal peak/standard deviation
                   sigma = 1/4.98,                          # susceptibility to primary infection 
                   d1 = 0.9325575,                          # relative susceptibility to secondary infection, relative to  primary infection
                   d2 = 0.7280313 ,                         # relative susceptibility to tertiary infection, relative to secondary infection
                   d3 = 0.35790461,                         # relative susceptibility to subsequent infections after third infection, relative to tertiary infection
                   Na = Na,                                 # population size for age group a
                   alpha = 0.5215113,                       # reduction in infectiousness of asymptomatic infections
                   p_contact = resceudata$contactMatrixPhy, # matrix total number of daily physical contacts made by age group a with age group b
                   c_contact = resceudata$contactMatrixCon, # matrix total number of daily conversational contacts made by age group a with age group b
                   propR = propR,                           # proportion of neonates born with protection
                   prob = prob,                             # proportion of persons in age group a who still have maternal protection
                   mu = 1863,                               # daily birth rate
                   xi = 1/133.5,                            # rate of loss maternal immunity
                   probA = probA,                           # probability infection is asymptomatic
                   g0 = 7.069876,                           # rate of loss of primary infectiousness
                   g1 = 0.9118119,                          # proportional decrease between secondary and primary infection
                   g2 = 0.6258454,                          # proportional decrease between tertiary and secondary infection
                   g3 = 1,                                  # duration of infectiousness for all infections after tertiary infection
                   omega = 1/358.9,                         # rate of loss of post-infection immunity
                   l1 = 0.08156498,                         # proportion of persons infected 
                   l2 = 0.5454382,                          # proportion of persons uninfected and protected
                   qp = 0.4904056,                          # probability of transmission per physical contact 
                   qc = 0.5083735                           # probability of transmission per conversational contact
)
############################ delta and gamma #####################
calculate_delta_gamma <- function(d1, d2, d3, g0, g1, g2) {
  delta1 <- d1
  delta2 <- d2*d1
  delta3 <- d1*d2*d3
  gamma0 <- 1/g0
  gamma1 <- 1/(g0*g1)
  gamma2 <- 1/(g0*g1*g2)
  
  return(list(delta1, delta2, delta3, gamma0, gamma1, gamma2))
} 
######################### initial proportion exposure ############## 
initialProportionExposure <- function(l,  #
                                      a1, # 
                                      a2  #
                                      ) {
  
  prop <- numeric(A)
  
  # proportions for the first three levels
  for (k in 0:2) {
    prop[k + 1] <- abs(ppois(k, l * a2) - ppois(k, l * a1)) / ((a2 - a1) * l)
  }
  
  # Fourth level
  prop[4] <- 1 - sum(prop[1:3])
  
  # vector of proportions
  return(prop)
}
############## initial M ############
initial_M <- function(parameter, A, ageStratification, populationPerAgeGroup) {
  xi <- 1/xi
  
  init_con <- numeric(NULL) 
  
  for (i in 1:(A - 1)) {
    ageStart <- ageStratification[i] * 365
    ageEnd <- ageStratification[i + 1] * 365
    
    cdfStart <- pexp(ageStart, rate = xi)
    cdfEnd <- pexp(ageEnd, rate = xi)
    
    init_con_temp <- (cdfEnd - cdfStart) / ((ageEnd - ageStart) * xi)
    init_con[i] <- init_con_temp * populationPerAgeGroup[i]
  }
  
  ageLastStart <- ageStratification[A] * 365
  ageLastEnd <- 90 * 365
  
  cdfLastStart <- pexp(ageLastStart, rate = xi)
  cdfLastEnd <- pexp(ageLastEnd, rate = xi)
  
  init_con_last <- (cdfLastEnd - cdfLastStart) / ((ageLastEnd - ageLastStart) * xi)
  init_con[A] <- init_con_last * populationPerAgeGroup[A]
  
  return(init_con)
}


################ initial states ################
generateInitialStates <- function(cov_c, A, populationPerAgeGroup, p_mat, pVHR, pHR, pLR, pA, parameters) {
  initialStates <- c()
  populationMatPro <- initial_M()  
  
  I1 <- I1
  I2 <- I2
  I3 <- 0.5
  
  xi <- 1.0 / xi
  gamma0 <- gamma0
  gamma1 <- gamma1
  gamma2 <- gamma2
  delta1 <- delta1
  delta2 <- delta2 
  delta3 <- delta3 
  
  for (a in 0:(A-1)) {
    if (a < A-1) {
      a1 <- ageStratification(a)  
      a2 <- ageStratification(a+1)  
    } else {
      a1 <- ageStratification(a)  
      a2 <- 90
    }
    
    prop <- initialProportionExposure(I3, a1, a2)
    pI1 <- prop[1]
    pI2 <- prop[2]
    pI3 <- prop[3]
    pI4 <- prop[4]
    
    age_size <- populationPerAgeGroup[a+1] - populationMatPro[a+1]
    
    for (s in 0:5) {
      if (a < 12) {
        s_prop <- switch(s+1,
                         `0` = 0.0,
                         `1` = cov_c,
                         `2` = 1-cov_c,
                         `3` = 0.0,
                         `4` = cov_c,
                         `5` = 1-cov_c)
      } else {
        s_prop <- switch(s+1,
                         `0` = p_mat[a+1] * (1-cov_c),
                         `1` = p_mat[a+1] * cov_c,
                         `2` = 1 - p_mat[a+1],
                         `3` = p_mat[a+1] * (1-cov_c),
                         `4` = p_mat[a+1] * cov_c,
                         `5` = 1 - p_mat[a+1])
      }
      
      for (r in 0:2) {
        r_prop <- switch(r+1,
                         `0` = pVHR[a+1],
                         `1` = pHR[a+1],
                         `2` = pLR[a+1])
##### my code #######      
        initialStates <- c(# maternal level 
                           N=8867713,
                           M=N*prob*xi,
                           
                           # exposure level 0
                           S0=N*(1-prob*xi)*prop[1]*(1-l1)*(1-l2), 
                           E0=N*(1-prob*xi)*prop[1]*(sigma/(gamma0+sigma))*l1, 
                           A0=N*(1-prob*xi)*prop[1]*(gamma0/(gamma0+sigma))*prob*l1,
                           I0=N*(1-prob*xi)*prop[1]*(gamma0/(gamma0+sigma))*(1-prob)*l1,
                           R0=N*(1-prob*xi)*prop[1]*(1-l1)*l2,
                           
                           # exposure level 1
                           S1=N*(1-prob*xi)*prop[2]*(1-(l1*delta1))*(1-l2), 
                           E1=N*(1-prob*xi)*prop[2]*(sigma/(gamma1+sigma))*l1*delta1, 
                           A1=N*(1-prob*xi)*prop[2]*(gamma1/(gamma1+sigma))*prob*l1*delta1,
                           I1=N*(1-prob*xi)*prop[2]*(gamma1/(gamma1+sigma))*(1-prob)*l1*delta1,
                           R1=N*(1-prob*xi)*prop[2]*(1-(l1*delta1))*l2,
                           
                           # exposure level 2
                           S2=N*(1-prob*xi)*prop[3]*(1-(l1*delta2))*(1-l2), 
                           E2=N*(1-prob*xi)*prop[3]*(sigma/(gamma2+sigma))*l1*delta2, 
                           A2=N*(1-prob*xi)*prop[3]*(gamma2/(gamma2+sigma))*prob*l1*delta2,
                           I2=N*(1-prob*xi)*prop[3]*(gamma2/(gamma2+sigma))*(1-prob)*l1*delta2,
                           R2=N*(1-prob*xi)*prop[3]*(1-(l1*delta2))*l2,
                           
                           # exposure level 3
                           S3=N*(1-prob*xi)*prop[4]*(1-(l1*delta3))*(1-l2), 
                           E3=N*(1-prob*xi)*prop[4]*(sigma/(gamma3+sigma))*l1*delta3, 
                           A3=N*(1-prob*xi)*prop[4]*(gamma3/(gamma3+sigma))*prob*l1*delta3,
                           I3=N*(1-prob*xi)*prop[4]*(gamma3/(gamma3+sigma))*(1-prob)*l1*delta3,
                           R3=N*(1-prob*xi)*prop[4]*(1-(l1*delta3))*l2,
                           
                           Z=0
        )
      }
    }

    initialStates <- c(initialStates, rep(0, 23))
  }
  
  return(initialStates)
}

