# Initialize the summation matrix
summation_matrix <- matrix(0, nrow = 5, ncol = 4)

# Initialize sums array to store N_b values
N_b <- numeric(5)

# Loop over each age group 'a' and infection level 'i'
for (i in 1:4) {
  for (a in 1:5) {
    # Initialize the summation for the current age group and infection level
    summation <- 0
    
    # Reset N_b calculation for each 'b' within the current loop
    for (b in 1:5) {
      # Calculate N_b for current 'b' and store in sums
      N_bi <- 0
      for (c in 0:3) {
        S_name <- paste0("S", c, b)
        matching_elements <- state_initial_age[grepl(S_name, names(state_initial_age))]
        if (length(matching_elements) > 0) {
          N_bi <- N_bi + sum(matching_elements)
        }
      }
      N_b[b] <- N_bi  # Store the computed N_b in sums array
      
      # Retrieve the A and I values for the current infection level 'i' and age group 'b'
      A_bi <- state_initial_age[paste0("A", i-1, b)]
      I_bi <- state_initial_age[paste0("I", i-1, b)]
      
      # Calculate the contribution of age group 'b' to the summation
      if (N_b[b] > 0) {  # Ensure that we do not divide by zero
        contribution <- ((contactMatrixPhy[a, b]) + (contactMatrixCon[a, b] * qc)) * (A_bi + I_bi) / N_b[b]
        summation <- summation + contribution
      }
    }
    
    # Store the result in the summation matrix
    summation_matrix[a, i] <- summation
  }
}

# The summation_matrix now contains the sum for each age group 'a' and infection level 'i'
summation_matrix
N_b
