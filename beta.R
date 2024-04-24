# define the age ranges
ranges <- list(c(1, 12), c(13, 16), c(17, 18), c(19, 23), c(24, 25))

# function to calculate average contact rate across ranges 
average_matrix_ranges <- function(data, ranges) {
  # initialize a matrix 
  averaged_matrix <- matrix(0, length(ranges), length(ranges))
  
  # averages for rows
  for (i in 1:length(ranges)) {
    row_range <- ranges[[i]]
    row_indices <- row_range[1]:row_range[2]
    
    #  averages for columns within each row block
    for (j in 1:length(ranges)) {
      col_range <- ranges[[j]]
      col_indices <- col_range[1]:col_range[2]
      
      # Extract the submatrix for the current row and column block
      submatrix <- data[row_indices, col_indices]
      
      # average for the submatrix
      averaged_matrix[i, j] <- mean(submatrix)
    }
  }
  
  return(averaged_matrix)
}

# get age-combined contact rates 
contactMatrixPhy <- average_matrix_ranges(uk_data_sum$contactMatrixPhy, ranges = ranges)
contactMatrixCon <- average_matrix_ranges(uk_data_sum$contactMatrixCon, ranges = ranges)


# initialize list to store lambda for each infection level and age group
lambda <- matrix(0, nrow = 4, ncol = 5)  
delta_i <- c(1,0.8908731,0.7912774,0.3428649)
qp = 0.09711828575 
qc = 0.99872580464


# Preallocate the lambda matrix
lambda_matrix <- matrix(0, nrow = 5, ncol = 4)

# Loop over each age group 'a' and infection level 'i' to calculate lambda
for (a in 1:5) {
  for (i in 1:4) {
    # Calculate the product of delta terms up to the infection level i
    delta_product <- prod(delta_i[1:i])
    
    # Calculate the time-varying factor
    beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))
    
    # Calculate lambda for the current age group 'a' and infection level 'i'
    lambda_matrix[a, i] <- qp * beta_value * delta_product * summation_matrix[a, i]
  }
}

# lambda_matrix now contains the lambda values for each age group 'a' and infection level 'i'
lambda_matrix
