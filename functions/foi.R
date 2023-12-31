############################### Age-dependent FOI ###################################
get_foi <- function(parameters, i, A, I, t, delta1, gamma0) {
  
  # define vectors and sequences to be used
  delta <- c(delta1, delta2, delta3)
  ages <- 1:25
  
  # output 
  lambda <- array(0, dim = c(t, ages, delta))  
  
  for (i in 1:length(delta)) {
    
    # seasonal forcing 
    lambda[a, i] <- qp * (1 + b1 * exp((t[d] - phi) / (2 * psi^2))) # beta seasonal variation 
    
    # iterate over age groups 
    sum_term_over_a <- 0 # sum of infected at t=0 for a=1
    
    # sum of infected 
    for (a in ages) {
      sum_term_over_a <- sum_term_over_a + (p_contact + (qc * c_contact)) / Na  * (alpha * A + I)
    }
    
    lambda[a, i] <- lambda[a, i]  * sum_term_over_a
    
  # iterate levels of delta over levels of exposure
  for (i_prime in 1:i) {
        
        lambda[a, i] <- lambda[a, i] * delta[i_prime]
      }
  }
  
  return(lambda)
}


