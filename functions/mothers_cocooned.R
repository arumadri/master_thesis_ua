####
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
