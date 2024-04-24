calculate_foi <- function(state_initial_age, contactMatrixPhy, contactMatrixCon, qp, b1, phi, psi, delta_i, qc) {
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
        
        # mixing 
        if (N_b[b] > 0) {  # avoid division by zero
          contacts <- ((contactMatrixPhy[a, b]) + (contactMatrixCon[a, b] * qc)) * (A_bi + I_bi) / N_b[b]
          summation <- summation + contacts
        }
      }
      
      contact_rate[a, i] <- summation
      
      # product of delta terms up to the infection level 'i'
      delta_product <- prod(delta_i[1:i])
      
      # Calculate the time-varying factor
      beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))
      
      # Calculate lambda for the current age group 'a' and infection level 'i'
      lambda_matrix[a, i] <- qp * beta_value * delta_product * summation
    }
  }
  
  return(lambda_matrix)
}


