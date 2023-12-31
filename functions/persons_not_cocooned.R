####
get_nwn <- function(prop_c, s, A, nwn_p, nwn_c) {
  nwn_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:A) {
      for (j in 1:A) {
        nwn_e[i, j] <- nwn_p(i, j)
      }
      for (j in 1:12) {
        nwn_e[i, j] <- nwn_p(i, j) * (1 - prop_c)
      }
    }
  } else if (s == 'c') {
    for (i in 1:A) {
      for (j in 1:A) {
        nwn_e[i, j] <- nwn_c(i, j)
      }
      for (j in 1:12) {
        nwn_e[i, j] <- nwn_c(i, j) * (1 - prop_c)
      }
    }
  } else {
    stop("Error cont")
  }
  
  return(nwn_e)
}
####
get_nwp <- function(prop_c, s, A, cnt_matrix_p, cnt_matrix_p_h, cnt_matrix_c, cnt_matrix_c_h, p_mat) {
  nwp_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:12) {
      for (j in 18:21) {
        nwp_e[i, j] <- (cnt_matrix_p_h[i, j] / 2.0) + (cnt_matrix_p[i, j] - cnt_matrix_p_h[i, j]) * p_mat[j] * (1 - prop_c)
      }
    }
    for (i in 13:25) {
      for (j in 18:21) {
        nwp_e[i, j] <- cnt_matrix_p[i, j] * p_mat[j] * (1 - prop_c)
      }
    }
  } else if (s == 'c') {
    for (i in 1:12) {
      for (j in 18:21) {
        nwp_e[i, j] <- (cnt_matrix_c_h[i, j] / 2.0) + (cnt_matrix_c[i, j] - cnt_matrix_c_h[i, j]) * p_mat[j] * (1 - prop_c)
      }
    }
    for (i in 13:25) {
      for (j in 18:21) {
        nwp_e[i, j] <- cnt_matrix_c[i, j] * p_mat[j] * (1 - prop_c)
      }
    }
  } else {
    stop("Error cont")
  }
  
  return(nwp_e)
}
####
get_nwc <- function(prop_c, s, A, nwn_p, cnt_matrix_p, cnt_matrix_p_h, p_mat, nwn_c, cnt_matrix_c, cnt_matrix_c_h) {
  nwc_e <- matrix(0, nrow = A, ncol = A)
  
  if (s == 'p') {
    for (i in 1:A) {
      for (j in 1:12) {
        nwc_e[i, j] <- nwn_p(i, j) * prop_c
      }
    }
    for (i in 1:12) {
      for (j in 18:21) {
        nwc_e[i, j] <- (cnt_matrix_p[i, j] - cnt_matrix_p_h[i, j]) * p_mat[j] * prop_c
      }
    }
    for (i in 13:25) {
      for (j in 18:21) {
        nwc_e[i, j] <- cnt_matrix_p[i, j] * p_mat[j] * prop_c
      }
    }
  } else if (s == 'c') {
    for (i in 1:A) {
      for (j in 1:12) {
        nwc_e[i, j] <- nwn_c(i, j) * prop_c
      }
    }
    for (i in 1:12) {
      for (j in 18:21) {
        nwc_e[i, j] <- (cnt_matrix_c[i, j] - cnt_matrix_c_h[i, j]) * p_mat[j] * prop_c
      }
    }
    for (i in 13:25) {
      for (j in 18:21) {
        nwc_e[i, j] <- cnt_matrix_c[i, j] * p_mat[j] * prop_c
      }
    }
  } else {
    stop("Error cont")
  }
  
  return(nwc_e)
}
