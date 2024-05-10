calculate_foi <- function(t1, S, A, I, delta_i, params) {

  lambda_matrix <- matrix(0, nrow = 5, ncol = 4)
  
  contactMatrixPhy <- params[['contactMatrixPhy']]
  contactMatrixCon <- params[['contactMatrixCon']]
  qp <- params[['qp']]
  qc <- params[['qc']]
  b1 <- params[['b1']]
  phi <- params[['phi']]
  psi <- params[['psi']]
  
  for (i in 1:4) {
    for (a in 1:5) {
      
      summation <- 0
      
      for (b in 1:5) {
        
        N_b = S[b+1]
        
        if (N_b > 0) {  # avoid division by zero
          contacts = ((contactMatrixPhy[a, b] + contactMatrixCon[a, b] * qc) * (A[b+1] + I[b+1])) / N_b
          summation = contacts
        }
      }
      
      # delta product up to infection level 'i'
      delta_product <- prod(delta_i[1:i])
      
      # beta
      beta_value <- 1 + b1 * exp(-((t1 - phi)^2) / (2 * psi^2))
      
      lambda_matrix[a, i] <- qp * beta_value * delta_product * summation
    }
  }
  
  return(lambda_matrix)
}
