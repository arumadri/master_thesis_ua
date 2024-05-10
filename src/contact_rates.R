# Initialize the summation matrix
summation_matrix <- matrix(0, nrow = 5, ncol = 4)

N_b <- numeric(5)

for (i in 1:4) {
  for (a in 1:5) {
    
    summation <- 0
    
    for (b in 1:5) {
      
      N_bi <- 0
      for (c in 0:3) {
        S_name <- paste0("S", c, b)
        matching_elements <- state_initial_age[grepl(S_name, names(state_initial_age))]
        if (length(matching_elements) > 0) {
          N_bi <- N_bi + sum(matching_elements)
        }
      }
      N_b[b] <- N_bi  
      
      A_bi <- state_initial_age[paste0("A", i-1, b)]
      I_bi <- state_initial_age[paste0("I", i-1, b)]
      
      if (N_b[b] > 0) { # do not divide by zero
        contribution <- ((contactMatrixPhy[a, b]) + (contactMatrixCon[a, b] * qc)) * (A_bi + I_bi) / N_b[b]
        summation <- summation + contribution
      }
    }
    
    summation_matrix[a, i] <- summation
  }
}

# The summation_matrix now contains the sum for each age group 'a' and infection level 'i'
summation_matrix
N_b
