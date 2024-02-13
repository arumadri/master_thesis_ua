##### ODE description 
ODE_desc <- list(
  # states and variables 
  M = 0, S0 = 0, S1 = 0, S2 = 0, S3 = 0, SN = 0, E0 = 0, E1 = 0, E2 = 0, E3 = 0, 
  A0 = 0, A1 = 0, A2 = 0, A3 = 0, I0 = 0, I1 = 0, I2 = 0, I3 = 0, R0 = 0, R1 = 0, R2 = 0, R3 = 0,
  M_ = 0, S0_ = 0, S1_ = 0, S2_ = 0, S3_ = 0, E0_ = 0, E1_ = 0, E2_ = 0, E3_ = 0,
  A0_ = 0, A1_ = 0, A2_ = 0, A3_ = 0, I0_ = 0, I1_ = 0, I2_ = 0, I3_ = 0, R0_ = 0, R1_ = 0, R2_ = 0, R3_ = 0, N = 0,
  Mv = 0, S0v = 0, S1v = 0, S2v = 0, S3v = 0, E0v = 0, E1v = 0, E2v = 0, E3v = 0, A0v = 0, A1v = 0, A2v = 0, A3v = 0, I0v = 0, I1v = 0, I2v = 0, I3v = 0, R0v = 0, R1v = 0, R2v = 0, R3v = 0,
  Mv_ = 0, S0v_ = 0, S1v_ = 0, S2v_ = 0, S3v_ = 0, E0v_ = 0, E1v_ = 0, E2v_ = 0, E3v_ = 0, A0v_ = 0, A1v_ = 0, A2v_ = 0, A3v_ = 0, I0v_ = 0, I1v_ = 0, I2v_ = 0, I3v_ = 0, R0v_ = 0, R1v_ = 0, R2v_ = 0, R3v_ = 0,
  xi = 0, si = 0, ga0 = 0, ga1 = 0, ga2 = 0, ga3 = 0, d1 = 0, d2 = 0, d3 = 0, a1 = 0, a2 = 0, a3 = 0, alpha_i = 0, rho = 0, om = 0, b1 = 0, qp = 0, qc = 0, psi = 0, phi = 0, beta = 0,
  
  A = 0,
  
  p_vul = 0,
  I_temp_p = 0, I_temp_c = 0, I_temp_n = 0, I_temp_p_v = 0, I_temp_c_v = 0, I_temp_n_v = 0,
  x_tot = 0, x_tot_1 = 0, x_tot_2 = 0, xi_b = 0,
  
  pVHR = numeric(), pHR = numeric(), pLR = numeric(), p_mat = numeric(),
  populationPerAgeGroup = numeric(), eta = numeric(), pA = numeric(),
  VHR_g = c(0.653012, 0.291842, 1, 1, 1, 0.0321774, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  LR_g = c(1.00126, 1.00168, 1, 1, 1, 1.00067, 1, 1, 1.00002, 1, 1, 1, 1, 1, 1, 1, 1),
  
  cal_pal = matrix(), cal_mAB_VHR = matrix(), cal_mAB_HR = matrix(), cal_mAB_LR = matrix(), 
  cal_LAV_HR = matrix(), cal_LAV_LR = matrix(), cal_mat_LR = matrix(),
  
  cal_pal_t = matrix(), cal_mAB_VHR_t = matrix(), cal_mAB_HR_t = matrix(), cal_mAB_LR_t = matrix(), 
  cal_LAV_HR_t = matrix(), cal_LAV_LR_t = matrix(), cal_mat_LR_t = matrix(),
  
  cal_pal_dose = matrix(), cal_mAB_VHR_dose = matrix(), cal_mAB_HR_dose = matrix(), cal_mAB_LR_dose = matrix(), 
  cal_LAV_HR_dose = matrix(), cal_LAV_LR_dose = matrix(), vac_cal_dose = matrix(),
  
  vac_calendar = list(), vac_info = list(), vac_dose = list(),
  
  cov_c = 0,  
  om_mab = 0,
  xi_boost = 0,
  
  direct = FALSE,
  
  u18p = 0, u19p = 0, u20p = 0,
  
  N_tot = numeric(A), N_tot_n = numeric(A), N_tot_c = numeric(A), N_tot_p = numeric(A),
  N_tot_n_v = numeric(A), N_tot_c_v = numeric(A), N_tot_p_v = numeric(A),
  
  prop_n = numeric(A), prop_c = numeric(A), prop_p = numeric(A),
  prop_nv = numeric(A), prop_cv = numeric(A), prop_pv = numeric(A),
  prop_empty = numeric(A),
  
  PS <- numeric()
)
#### ODE function 
ODE_desc <- function(finELL_t, vac_calendar_t, vac_dose_t, vac_info_t, cov_c_t) {
  
  source("functions/RunInterventions.R")
  
  finELL <- finELL_t
  vac_calendar <- vac_calendar_t
  vac_dose <- vac_dose_t
  vac_info <- vac_info_t
  cov_c <- cov_c_t
  
  # Perform calculations and assign additional properties
  xi <- 1.0 / finELL$parameterValues[["xi"]]
  si <- 1.0 / finELL$parameterValues[["si"]]
  ga0 <- 1.0 / finELL$parameterValues[["ga0"]]
  ga1 <- 1.0 / (finELL$parameterValues[["ga0"]] * finELL$parameterValues[["g1"]])
  ga2 <- 1.0 / (finELL$parameterValues[["ga0"]] * finELL$parameterValues[["g1"]] * finELL$parameterValues[["g2"]])
  ga3 <- ga2
  om <- 1.0 / finELL$parameterValues[["om"]]
  
  rho <- 1.0
  alpha_i <- finELL$parameterValues[["alpha_i"]]
  d1 <- finELL$parameterValues[["d1"]]
  d2 <- finELL$parameterValues[["d1"]] * finELL$parameterValues[["d2"]]
  d3 <- finELL$parameterValues[["d1"]] * finELL$parameterValues[["d2"]] * finELL$parameterValues[["d3"]]
  a1 <- 1.0
  a2 <- 1.0
  a3 <- 1.0
  
  phi <- finELL$parameterValues[["phi"]]
  qp <- finELL$parameterValues[["qp"]]
  qc <- finELL$parameterValues[["qc"]]
  b1 <- finELL$parameterValues[["b1"]]
  psi <- finELL$parameterValues[["psi"]]
  
  A <- finELL$A
  eta <- finELL$eta
  
  u18p <- finELL$u_p[1]
  u19p <- finELL$u_p[2]
  u20p <- finELL$u_p[3]
  
  populationPerAgeGroup <- finELL$populationPerAgeGroup
  pVHR <- finELL$pVHR
  pHR <- finELL$pHR
  pLR <- finELL$pLR
  p_mat <- finELL$p_mat

  dailyBirthRate <- finELL$dailyBirthRate
  pA <- finELL$pA
  
  I_n <- rep(0, A)
  I_c <- rep(0, A)
  I_p <- rep(0, A)
  I_n_v <- rep(0, A)
  I_c_v <- rep(0, A)
  I_p_v <- rep(0, A)
  
  N_tot <- rep(0, A)
  N_tot_n <- rep(0, A)
  N_tot_c <- rep(0, A)
  N_tot_p <- rep(0, A)
  N_tot_n_v <- rep(0, A)
  N_tot_c_v <- rep(0, A)
  N_tot_p_v <- rep(0, A)
  
  prop_n <- rep(0, A)
  prop_c <- rep(0, A)
  prop_p <- rep(0, A)
  prop_nv <- rep(0, A)
  prop_cv <- rep(0, A)
  prop_pv <- rep(0, A)
  prop_empty <- rep(0, A)
  
  PS <- rep(0, (A - 1))
  
  cal_pal_t <- vac_calendar[["pal"]]
  cal_mAB_VHR_t <- vac_calendar[["mAB_VHR"]]
  cal_mAB_HR_t <- vac_calendar[["mAB_HR"]]
  cal_mAB_LR_t <- vac_calendar[["mAB_LR"]]
  cal_LAV_HR_t <- vac_calendar[["LAV_HR"]]
  cal_LAV_LR_t <- vac_calendar[["LAV_LR"]]
  cal_mat_LR_t <- vac_calendar[["mat_LR"]]
  
  cal_pal_dose <- vac_dose[["pal"]]
  cal_mAB_VHR_dose <- vac_dose[["mAB_VHR"]]
  cal_mAB_HR_dose <- vac_dose[["mAB_HR"]]
  cal_mAB_LR_dose <- vac_dose[["mAB_LR"]]
  cal_LAV_HR_dose <- vac_dose[["LAV_HR"]]
  cal_LAV_LR_dose <- vac_dose[["LAV_LR"]]
  vac_cal_dose <- vac_dose[["mat_LR"]]
  
  direct < vac_info[["direct"]]
  om_mab <- vac_info[["om_mab"]]
  xi_boost <- vac_info[["xi_boost"]]

######## operator    
operator <- function(x, dxdt, t) {
    ag <- 455
    sg <- 72
    rg <- 24
    # CALENDARS
    vac_cal_vhr <- matrix(0, 365, A)  
    vac_cal <- matrix(0, 365, A)
    
    if (t < 365) {
      cal_pal <- matrix(0, 365, A)
      vac_cal <- matrix(0, 365, A)
      cal_mAB_VHR <- matrix(0, 365, A)
      cal_mAB_HR <- matrix(0, 365, A)
      cal_mAB_LR <- matrix(0, 365, A)
      cal_LAV_HR <- matrix(0, 365, A)
      cal_LAV_LR <- matrix(0, 365, A)
    } else {
      cal_pal <- cal_pal_t
      vac_cal <- cal_mat_LR_t
      cal_mAB_VHR <- cal_mAB_VHR_t
      cal_mAB_HR <- cal_mAB_HR_t
      cal_mAB_LR <- cal_mAB_LR_t
      cal_LAV_HR <- cal_LAV_HR_t
      cal_LAV_LR <- cal_LAV_LR_t
    }

##### populations
tot_coc <- 0
tot_coc_v <- 0
tot_inf <- 0
tot_inf_v <- 0

for (a in 1:A) {
  tot_temp_n <- 0
  tot_temp_c <- 0
  tot_temp_p <- 0
  tot_temp_nv <- 0
  tot_temp_cv <- 0
  tot_temp_pv <- 0
  tot_temp_vhr <- 0
  tot_temp_vhrv <- 0
  tot_temp_hr <- 0
  tot_temp_hrv <- 0
  tot_temp_lr <- 0
  tot_temp_lrv <- 0
  
  for (j in 1:24) {
    for (r in 1:3) {
      tot_temp_p <- tot_temp_p + x[(a-1)*ag + 0*sg + (r-1)*rg + j]
      tot_temp_c <- tot_temp_c + x[(a-1)*ag + 1*sg + (r-1)*rg + j]
      tot_temp_n <- tot_temp_n + x[(a-1)*ag + 2*sg + (r-1)*rg + j]
      tot_temp_pv <- tot_temp_pv + x[(a-1)*ag + 3*sg + (r-1)*rg + j]
      tot_temp_cv <- tot_temp_cv + x[(a-1)*ag + 4*sg + (r-1)*rg + j]
      tot_temp_nv <- tot_temp_nv + x[(a-1)*ag + 5*sg + (r-1)*rg + j]
    }
  }
  
  for (j in 1:24) {
    for (s in 1:3) {
      tot_temp_vhr <- tot_temp_vhr + x[(a-1)*ag + s*sg + 0*rg + j]
      tot_temp_hr <- tot_temp_hr + x[(a-1)*ag + s*sg + 1*rg + j]
      tot_temp_lr <- tot_temp_lr + x[(a-1)*ag + s*sg + 2*rg + j]
      tot_temp_vhrv <- tot_temp_vhrv + x[(a-1)*ag + (s+3)*sg + 0*rg + j]
      tot_temp_hrv <- tot_temp_hrv + x[(a-1)*ag + (s+3)*sg + 1*rg + j]
      tot_temp_lrv <- tot_temp_lrv + x[(a-1)*ag + (s+3)*sg + 2*rg + j]
    }
  }
  
  if (a <= 12) {
    tot_coc <- tot_coc + tot_temp_c
    tot_coc_v <- tot_coc_v + tot_temp_cv
    tot_inf <- tot_inf + tot_temp_c + tot_temp_n
    tot_inf_v <- tot_inf_v + tot_temp_cv + tot_temp_nv
  }
  
  N_tot_n[a] <- tot_temp_n
  N_tot_c[a] <- tot_temp_c
  N_tot_p[a] <- tot_temp_p
  N_tot_n_v[a] <- tot_temp_nv
  N_tot_c_v[a] <- tot_temp_cv
  N_tot_p_v[a] <- tot_temp_pv
}

### force of infection
# Define phi_c and phi_c_v
phi_c <- cov_c
phi_c_v <- cov_c

# 2.1 Get contact matrices
cnt_matrix_cwc_p <- get_cwc(phi_c, 'p')
cnt_matrix_cwp_p <- get_cwp(phi_c, 'p')
cnt_matrix_cwn_p <- get_cwn(phi_c, 'p')
cnt_matrix_pwc_p <- get_pwc(phi_c, 'p')
cnt_matrix_pwp_p <- get_pwp(phi_c, 'p')
cnt_matrix_pwn_p <- get_pwn(phi_c, 'p')
cnt_matrix_nwc_p <- get_nwc(phi_c, 'p')
cnt_matrix_nwp_p <- get_nwp(phi_c, 'p')
cnt_matrix_nwn_p <- get_nwn(phi_c, 'p')

cnt_matrix_cwc_c <- get_cwc(phi_c, 'c')
cnt_matrix_cwp_c <- get_cwp(phi_c, 'c')
cnt_matrix_cwn_c <- get_cwn(phi_c, 'c')
cnt_matrix_pwc_c <- get_pwc(phi_c, 'c')
cnt_matrix_pwp_c <- get_pwp(phi_c, 'c')
cnt_matrix_pwn_c <- get_pwn(phi_c, 'c')
cnt_matrix_nwc_c <- get_nwc(phi_c, 'c')
cnt_matrix_nwp_c <- get_nwp(phi_c, 'c')
cnt_matrix_nwn_c <- get_nwn(phi_c, 'c')

cnt_matrix_cwc_p_v <- get_cwc(phi_c_v, 'p')
cnt_matrix_cwp_p_v <- get_cwp(phi_c_v, 'p')
cnt_matrix_cwn_p_v <- get_cwn(phi_c_v, 'p')
cnt_matrix_pwc_p_v <- get_pwc(phi_c_v, 'p')
cnt_matrix_pwp_p_v <- get_pwp(phi_c_v, 'p')
cnt_matrix_pwn_p_v <- get_pwn(phi_c_v, 'p')
cnt_matrix_nwc_p_v <- get_nwc(phi_c_v, 'p')
cnt_matrix_nwp_p_v <- get_nwp(phi_c_v, 'p')
cnt_matrix_nwn_p_v <- get_nwn(phi_c_v, 'p')

cnt_matrix_cwc_c_v <- get_cwc(phi_c_v, 'c')
cnt_matrix_cwp_c_v <- get_cwp(phi_c_v, 'c')
cnt_matrix_cwn_c_v <- get_cwn(phi_c_v, 'c')
cnt_matrix_pwc_c_v <- get_pwc(phi_c_v, 'c')
cnt_matrix_pwp_c_v <- get_pwp(phi_c_v, 'c')
cnt_matrix_pwn_c_v <- get_pwn(phi_c_v, 'c')
cnt_matrix_nwc_c_v <- get_nwc(phi_c_v, 'c')
cnt_matrix_nwp_c_v <- get_nwp(phi_c_v, 'c')
cnt_matrix_nwn_c_v <- get_nwn(phi_c_v, 'c')

# 2.2 Seasonal forcing
t1 <- (as.integer(t) %% 365)
beta <- (1 + b1 * (1 + exp(-((t1 / 365.0 - phi)) * ((t1 / 365.0 - phi)) / (2 * psi * psi))))

# 3. DYNAMIC MATERNAL PROTECTION
# num_vec pos_wcb_ag = c(18, 21)
sum_wcb <- 0
sum_wcb_v <- 0
CB2_temp <- CB2_temp_v <- 0

for (a in 18:20) {
  for (r in 1:3) {
    CB2_temp <- CB2_temp + (x[a * ag + 2 * sg + r * rg + 1] +
                              x[a * ag + 2 * sg + r * rg + 6] +
                              x[a * ag + 2 * sg + r * rg + 11] +
                              x[a * ag + 2 * sg + r * rg + 16] +
                              x[a * ag + 2 * sg + r * rg + 2] +
                              x[a * ag + 2 * sg + r * rg + 7] +
                              x[a * ag + 2 * sg + r * rg + 12] +
                              x[a * ag + 2 * sg + r * rg + 17])
    
    CB2_temp_v <- CB2_temp_v + (x[a * ag + 5 * sg + r * rg + 1] +
                                  x[a * ag + 5 * sg + r * rg + 6] +
                                  x[a * ag + 5 * sg + r * rg + 11] +
                                  x[a * ag + 5 * sg + r * rg + 16] +
                                  x[a * ag + 5 * sg + r * rg + 2] +
                                  x[a * ag + 5 * sg + r * rg + 7] +
                                  x[a * ag + 5 * sg + r * rg + 12] +
                                  x[a * ag + 5 * sg + r * rg + 17])
  }
}

sum_wcb <- CB2_temp / sum(N_tot_n[18:20])
sum_wcb_v <- CB2_temp_v / sum(N_tot_n[18:20])

tot_C <- sum(N_tot_c[18:20])
tot_C_inv <- ifelse(tot_C < 1.0, 0, 1.0 / tot_C)

### ODES

# Initialization
protectpal <- 0
protectmabs <- 0
protectLAV <- 0
protectmat <- 0

# ODEs
for (a in 1:A) {
  r_prop <- 0
  protectpal <- 0
  protectmabs <- 0
  protectLAV <- 0
  protectmat <- 0
  I_temp_c <- 0
  I_temp_p <- 0
  I_temp_n <- 0
  
  for (k in 1:A) {
    if (N_tot_n[k] < 0.1) {
      N_tot_n_inv <- 0
    } else {
      N_tot_n_inv <- 1.0 / N_tot_n[k]
    }
    if (N_tot_c[k] < 0.1) {
      N_tot_c_inv <- 0
    } else {
      N_tot_c_inv <- 1.0 / N_tot_c[k]
    }
    if (N_tot_p[k] < 0.1) {
      N_tot_p_inv <- 0
    } else {
      N_tot_p_inv <- 1.0 / N_tot_p[k]
    }
    
    for (r in 1:3) {
      I_temp_p <- I_temp_p + (x[k * ag + 0 * sg + r * rg + 3] * alpha_i + x[k * ag + 0 * sg + r * rg + 4] +
                                a1 * (x[k * ag + 0 * sg + r * rg + 8] * alpha_i + x[k * ag + 0 * sg + r * rg + 9]) +
                                a2 * (x[k * ag + 0 * sg + r * rg + 13] * alpha_i + x[k * ag + 0 * sg + r * rg + 14]) +
                                a3 * (x[k * ag + 0 * sg + r * rg + 18] * alpha_i + x[k * ag + 0 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_pwp_p[a][k] + qc * cnt_matrix_pwp_c[a][k])) * N_tot_p_inv +
        (x[k * ag + 1 * sg + r * rg + 3] * alpha_i + x[k * ag + 1 * sg + r * rg + 4] +
           a1 * (x[k * ag + 1 * sg + r * rg + 8] * alpha_i + x[k * ag + 1 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 1 * sg + r * rg + 13] * alpha_i + x[k * ag + 1 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 1 * sg + r * rg + 18] * alpha_i + x[k * ag + 1 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_pwc_p[a][k] + qc * cnt_matrix_pwc_c[a][k])) * N_tot_c_inv +
        (x[k * ag + 2 * sg + r * rg + 3] * alpha_i + x[k * ag + 2 * sg + r * rg + 4] +
           a1 * (x[k * ag + 2 * sg + r * rg + 8] * alpha_i + x[k * ag + 2 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 2 * sg + r * rg + 13] * alpha_i + x[k * ag + 2 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 2 * sg + r * rg + 18] * alpha_i + x[k * ag + 2 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_pwn_p[a][k] + qc * cnt_matrix_pwn_c[a][k])) * N_tot_n_inv
      
      I_temp_c <- I_temp_c + (x[k * ag + 0 * sg + r * rg + 3] * alpha_i + x[k * ag + 0 * sg + r * rg + 4] +
                                a1 * (x[k * ag + 0 * sg + r * rg + 8] * alpha_i + x[k * ag + 0 * sg + r * rg + 9]) +
                                a2 * (x[k * ag + 0 * sg + r * rg + 13] * alpha_i + x[k * ag + 0 * sg + r * rg + 14]) +
                                a3 * (x[k * ag + 0 * sg + r * rg + 18] * alpha_i + x[k * ag + 0 * sg + r * rg + 19])) *
        (qc * (cnt_matrix_cwp_p[a][k] + qp * cnt_matrix_cwp_c[a][k])) * N_tot_p_inv +
        (x[k * ag + 1 * sg + r * rg + 3] * alpha_i + x[k * ag + 1 * sg + r * rg + 4] +
           a1 * (x[k * ag + 1 * sg + r * rg + 8] * alpha_i + x[k * ag + 1 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 1 * sg + r * rg + 13] * alpha_i + x[k * ag + 1 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 1 * sg + r * rg + 18] * alpha_i + x[k * ag + 1 * sg + r * rg + 19])) *
        (qc * (cnt_matrix_cwc_p[a][k] + qp * cnt_matrix_cwc_c[a][k])) * N_tot_c_inv +
        (x[k * ag + 2 * sg + r * rg + 3] * alpha_i + x[k * ag + 2 * sg + r * rg + 4] +
           a1 * (x[k * ag + 2 * sg + r * rg + 8] * alpha_i + x[k * ag + 2 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 2 * sg + r * rg + 13] * alpha_i + x[k * ag + 2 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 2 * sg + r * rg + 18] * alpha_i + x[k * ag + 2 * sg + r * rg + 19])) *
        (qc * (cnt_matrix_cwn_p[a][k] + qp * cnt_matrix_cwn_c[a][k])) * N_tot_n_inv
      
      I_temp_n <- I_temp_n + (x[k * ag + 0 * sg + r * rg + 3] * alpha_i + x[k * ag + 0 * sg + r * rg + 4] +
                                a1 * (x[k * ag + 0 * sg + r * rg + 8] * alpha_i + x[k * ag + 0 * sg + r * rg + 9]) +
                                a2 * (x[k * ag + 0 * sg + r * rg + 13] * alpha_i + x[k * ag + 0 * sg + r * rg + 14]) +
                                a3 * (x[k * ag + 0 * sg + r * rg + 18] * alpha_i + x[k * ag + 0 * sg + r * rg + 19])) *
        (qn * (cnt_matrix_nwp_p[a][k] + qc * cnt_matrix_nwp_c[a][k])) * N_tot_p_inv +
        (x[k * ag + 1 * sg + r * rg + 3] * alpha_i + x[k * ag + 1 * sg + r * rg + 4] +
           a1 * (x[k * ag + 1 * sg + r * rg + 8] * alpha_i + x[k * ag + 1 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 1 * sg + r * rg + 13] * alpha_i + x[k * ag + 1 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 1 * sg + r * rg + 18] * alpha_i + x[k * ag + 1 * sg + r * rg + 19])) *
        (qn * (cnt_matrix_nwc_p[a][k] + qc * cnt_matrix_nwc_c[a][k])) * N_tot_c_inv +
        (x[k * ag + 2 * sg + r * rg + 3] * alpha_i + x[k * ag + 2 * sg + r * rg + 4] +
           a1 * (x[k * ag + 2 * sg + r * rg + 8] * alpha_i + x[k * ag + 2 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 2 * sg + r * rg + 13] * alpha_i + x[k * ag + 2 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 2 * sg + r * rg + 18] * alpha_i + x[k * ag + 2 * sg + r * rg + 19])) *
        (qn * (cnt_matrix_nwn_p[a][k] + qc * cnt_matrix_nwn_c[a][k])) * N_tot_n_inv
    }
  }

  I_temp_c_v <- 0.0
  I_temp_p_v <- 0.0
  I_temp_n_v <- 0.0
  
  for (k in 1:A) {
    N_tot_n_v_inv <- ifelse(N_tot_n_v[k] < 1, 0, 1.0 / N_tot_n_v[k])
    N_tot_c_v_inv <- ifelse(N_tot_c_v[k] < 1, 0, 1.0 / N_tot_c_v[k])
    N_tot_p_v_inv <- ifelse(N_tot_p_v[k] < 1, 0, 1.0 / N_tot_p_v[k])
    
    for (r in 1:3) {
      I_temp_p_v <- I_temp_p_v + (x[k * ag + 3 * sg + r * rg + 3] * alpha_i +
                                    x[k * ag + 3 * sg + r * rg + 4] + a1 * (x[k * ag + 3 * sg + r * rg + 8] * alpha_i + x[k * ag + 3 * sg + r * rg + 9]) +
                                    a2 * (x[k * ag + 3 * sg + r * rg + 13] * alpha_i + x[k * ag + 3 * sg + r * rg + 14]) +
                                    a3 * (x[k * ag + 3 * sg + r * rg + 18] * alpha_i + x[k * ag + 3 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_pwp_p_v[a][k] + qc * cnt_matrix_pwp_c_v[a][k])) * N_tot_p_v_inv +
        (x[k * ag + 4 * sg + r * rg + 3] * alpha_i +
           x[k * ag + 4 * sg + r * rg + 4] + a1 * (x[k * ag + 4 * sg + r * rg + 8] * alpha_i + x[k * ag + 4 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 4 * sg + r * rg + 13] * alpha_i + x[k * ag + 4 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 4 * sg + r * rg + 18] * alpha_i + x[k * ag + 4 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_pwc_p_v[a][k] + qc * cnt_matrix_pwc_c_v[a][k])) * N_tot_c_v_inv +
        (x[k * ag + 5 * sg + r * rg + 3] * alpha_i +
           x[k * ag + 5 * sg + r * rg + 4] + a1 * (x[k * ag + 5 * sg + r * rg + 8] * alpha_i + x[k * ag + 5 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 5 * sg + r * rg + 13] * alpha_i + x[k * ag + 5 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 5 * sg + r * rg + 18] * alpha_i + x[k * ag + 5 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_pwn_p_v[a][k] + qc * cnt_matrix_pwn_c_v[a][k])) * N_tot_n_v_inv;
      
      I_temp_c_v <- I_temp_c_v + (x[k * ag + 3 * sg + r * rg + 3] * alpha_i +
                                    x[k * ag + 3 * sg + r * rg + 4] + a1 * (x[k * ag + 3 * sg + r * rg + 8] * alpha_i + x[k * ag + 3 * sg + r * rg + 9]) +
                                    a2 * (x[k * ag + 3 * sg + r * rg + 13] * alpha_i + x[k * ag + 3 * sg + r * rg + 14]) +
                                    a3 * (x[k * ag + 3 * sg + r * rg + 18] * alpha_i + x[k * ag + 3 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_cwp_p_v[a][k] + qc * cnt_matrix_cwp_c_v[a][k])) * N_tot_p_v_inv +
        (x[k * ag + 4 * sg + r * rg + 3] * alpha_i +
           x[k * ag + 4 * sg + r * rg + 4] + a1 * (x[k * ag + 4 * sg + r * rg + 8] * alpha_i + x[k * ag + 4 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 4 * sg + r * rg + 13] * alpha_i + x[k * ag + 4 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 4 * sg + r * rg + 18] * alpha_i + x[k * ag + 4 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_cwc_p_v[a][k] + qc * cnt_matrix_cwc_c_v[a][k])) * N_tot_c_v_inv +
        (x[k * ag + 5 * sg + r * rg + 3] * alpha_i +
           x[k * ag + 5 * sg + r * rg + 4] + a1 * (x[k * ag + 5 * sg + r * rg + 8] * alpha_i + x[k * ag + 5 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 5 * sg + r * rg + 13] * alpha_i + x[k * ag + 5 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 5 * sg + r * rg + 18] * alpha_i + x[k * ag + 5 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_cwn_p_v[a][k] + qc * cnt_matrix_cwn_c_v[a][k])) * N_tot_n_v_inv;
      
      I_temp_n_v <- I_temp_n_v + (x[k * ag + 3 * sg + r * rg + 3] * alpha_i +
                                    x[k * ag + 3 * sg + r * rg + 4] +
                                    a1 * (x[k * ag + 3 * sg + r * rg + 8] * alpha_i + x[k * ag + 3 * sg + r * rg + 9]) +
                                    a2 * (x[k * ag + 3 * sg + r * rg + 13] * alpha_i + x[k * ag + 3 * sg + r * rg + 14]) +
                                    a3 * (x[k * ag + 3 * sg + r * rg + 18] * alpha_i + x[k * ag + 3 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_nwp_p_v[a][k] + qc * cnt_matrix_nwp_c_v[a][k])) * N_tot_p_v_inv +
        (x[k * ag + 4 * sg + r * rg + 3] * alpha_i +
            x[k * ag + 4 * sg + r * rg + 4] +
            a1 * (x[k * ag + 4 * sg + r * rg + 8] * alpha_i + x[k * ag + 4 * sg + r * rg + 9]) +
            a2 * (x[k * ag + 4 * sg + r * rg + 13] * alpha_i + x[k * ag + 4 * sg + r * rg + 14]) +
            a3 * (x[k * ag + 4 * sg + r * rg + 18] * alpha_i + x[k * ag + 4 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_nwc_p_v[a][k] + qc * cnt_matrix_nwc_c_v[a][k])) * N_tot_c_v_inv +
        (x[k * ag + 5 * sg + r * rg + 3] * alpha_i +
           x[k * ag + 5 * sg + r * rg + 4] +
           a1 * (x[k * ag + 5 * sg + r * rg + 8] * alpha_i + x[k * ag + 5 * sg + r * rg + 9]) +
           a2 * (x[k * ag + 5 * sg + r * rg + 13] * alpha_i + x[k * ag + 5 * sg + r * rg + 14]) +
           a3 * (x[k * ag + 5 * sg + r * rg + 18] * alpha_i + x[k * ag + 5 * sg + r * rg + 19])) *
        (qp * (cnt_matrix_nwn_p_v[a][k] + qc * cnt_matrix_nwn_c_v[a][k])) * N_tot_n_v_inv
               
    }
  }
  
  if (direct) {
    I_temp_p_v <- I_temp_p
    I_temp_c_v <- I_temp_c
    I_temp_n_v <- I_temp_n
  }

  pj <- max(ag * (a - 1), 0)
  cj <- max(ag * a, 0)
  muBp <- 0
  muBc <- 0
  muBn <- 0
  muBpv <- 0
  muBcv <- 0
  muBnv <- 0
  mu_mat <- 0
  cl <- 0
  Icp <- 0
  Icc <- 0
  up <- 0
  kpc <- 0
  kpd <- 0
  pro <- 0
  pro_v_in <- 0
  pro_v_out <- 0
  rp <- 0
  p_vulp <- 0
  p_vulc <- 0
  p_vuln <- 0
  p_vulpv <- 0
  p_vulcv <- 0
  p_vulnv <- 0
  xi_bp <- 0
  xi_bc <- 0
  xi_bn <- 0
  xi_bpv <- 0
  xi_bcv <- 0
  xi_bnv <- 0
  ej1 <- eta[a + 1]
  ej <- eta[a]
  lossP <- 0
  lossMS0 <- 0
  lossMS1 <- 0
  
  # Birth rate into each social group
  if (a == 0) {
    p_vulp <- 0
    p_vulc <- sum_wcb
    p_vuln <- sum_wcb
    p_vulpv <- 0
    p_vulcv <- sum_wcb_v
    p_vulnv <- sum_wcb_v
    
    muBp <- 0
    muBc <- dailyBirthRate * (cov_c)
    muBn <- dailyBirthRate * (1 - cov_c)
    muBpv <- 0
    muBcv <- dailyBirthRate * (cov_c)
    muBnv <- dailyBirthRate * (1 - cov_c)
    
    xi_bp <- 1
    xi_bc <- 1
    xi_bn <- 1
    xi_bpv <- 1
    xi_bcv <- xi_boost
    xi_bnv <- 1
    
  } else {
    
    p_vulp <- 0
    p_vulc <- 0
    p_vuln <- 0
    p_vulpv <- 0
    p_vulcv <- 0
    p_vulnv <- 0
    
    muBp <- 0
    muBc <- 0
    muBn <- 0
    muBpv <- 0
    muBcv <- 0
    muBnv <- 0
    
    xi_bp <- 1
    xi_bc <- 1
    xi_bn <- 1
    xi_bpv <- 1
    xi_bcv <- 1
    xi_bnv <- 1
  }
  
  for (s in 1:6) {
    # Initializing variables
    kpd <- 0
    kpc <- 0
    pro <- 0
    pro_v_in <- 0
    pro_v_out <- 0
    cl <- 0
    Icp <- 1
    Icc <- 0
    up <- 0
    
    # Initializing additional variables
    u <- 0
    In <- 0
    mu <- 0
    vac_c_o <- 0
    vac_c_i <- 0
    kp <- 0
    
    if (a < 12) {
      if (s == 1) {
        kp <- 0 * sg
        u <- 0
        In <- I_temp_p
        mu <- muBp
        mu_mat <- 0 
        p_vul <- p_vulp 
        xi_b <-  xi_bp
      } else if (s == 2) {
        kp <- 1 * sg
        u <- cov_c
        In <- I_temp_c
        mu <- muBc
        mu_mat <- 0 
        p_vul <- p_vulc 
        xi_b <-  xi_bc
      } else if (s == 3) {
        kp <- 2 * sg
        u <- 1 - cov_c
        In <- I_temp_n
        mu <- muBn
        mu_mat <- 0 
        p_vul <- p_vuln 
        xi_b <-  xi_bn
      } else if (s == 4) {
        kp <- 3 * sg
        u <- 0
        In <- I_temp_p_v
        mu <- muBpv
        mu_mat <- 0 
        p_vul <- p_vulpv 
        xi_b <-  xi_bpv
      } else if (s == 5) {
        kp <- 4 * sg
        u <- cov_c
        In <- I_temp_c_v
        mu <- muBcv * (1 - vac_cal(t1, 0))
        mu_mat <- muBcv*vac_cal(t1,0)
        p_vul <- p_vulcv 
        xi_b <-  xi_bcv
      } else {
        kp <- 5 * sg
        u <- 1 - cov_c
        In <- I_temp_n_v
        mu <- muBnv
        mu_mat <- 0 
        p_vul <- p_vulnv
        xi_b <-  xi_bnv
      }
    } else {
      if (s == 1) {
        kp <- 0 * sg
        u <- p_mat[a] * (1 - cov_c)
        In <- I_temp_p
        mu <- 0
        mu_mat <- 0
        p_vul <- p_vulp
      } else if (s == 2) {
        kp <- 1 * sg
        u <- p_mat[a] * cov_c
        In <- I_temp_c
        mu <- 0
        mu_mat <- 0
        p_vul <- p_vulc
      } else if (s == 3) {
        kp <- 2 * sg
        u <- 1 - p_mat[a]
        In <- I_temp_n
        mu <- 0
        mu_mat <- 0
        p_vul <- p_vuln
      } else if (s == 4) {
        kp <- 3 * sg
        u <- p_mat[a] * (1 - cov_c)
        In <- I_temp_p_v
        mu <- 0
        mu_mat <- 0
        p_vul <- p_vulpv
      } else if (s == 5) {
        kp <- 4 * sg
        u <- p_mat[a] * cov_c
        In <- I_temp_c_v
        mu <- 0
        mu_mat <- 0
        p_vul <- p_vulcv
      } else {
        kp <- 5 * sg
        u <- 1 - p_mat[a]
        In <- I_temp_n_v
        mu <- 0
        mu_mat <- 0
        p_vul <- p_vulnv
      }
    }

    for (r in 1:3) {
      # Rcpp::Rcout << "Risk groups parameters: " 
      
      PST <- 0
      x_tot <- x_tot_1 <- x_tot_2 <- 0
      if (r == 1) {
        rp <- pVHR[a]
      } else if (r == 2) {
        rp <- pHR[a]
      } else if (r == 3) {
        rp <- pLR[a]
      } else {
        cat("ERROR\n")
      }
      
      if (r == 3) {
        # Rcpp::Rcout << "Place where r == 3." 
        if (a %in% c(18, 19, 20)) {
          up <- switch(u18p, u19p, u20p)
        }
      }
      
      if (s < 3) {
        # Rcpp::Rcout << "Place where s < 3." 
        if (a < 12) {
          for (i in 1:21) {
            PS[i] <- sum(x[pj + 0*sg + 0*rg + i] + x[pj + 0*sg + 1*rg + i] + x[pj + 0*sg + 2*rg + i] +
                           x[pj + 1*sg + 0*rg + i] + x[pj + 1*sg + 1*rg + i] + x[pj + 1*sg + 2*rg + i] +
                           x[pj + 2*sg + 0*rg + i] + x[pj + 2*sg + 1*rg + i] + x[pj + 2*sg + 2*rg + i])
          }
          
          for (i in 22:24) {
            PS[i] <- sum(x[pj + 0*sg + r*rg + i] + x[pj + 1*sg + r*rg + i] + x[pj + 2*sg + r*rg + i])
          }
        } else {
          for (i in 1:24) {
            PS[i] <- sum(x[pj + 0*sg + 0*rg + i] + x[pj + 0*sg + 1*rg + i] + x[pj + 0*sg + 2*rg + i] +
                           x[pj + 1*sg + 0*rg + i] + x[pj + 1*sg + 1*rg + i] + x[pj + 1*sg + 2*rg + i] +
                           x[pj + 2*sg + 0*rg + i] + x[pj + 2*sg + 1*rg + i] + x[pj + 2*sg + 2*rg + i])
          }
        }
      } else {
        # Rcpp::Rcout << "Place where s >= 3." 
        if (a < 12) {
          # Rcpp::Rcout << "Place where a < 12." 
          for (i in 1:21) {
            PS[i] <- sum(x[pj + 3*sg + 0*rg + i] + x[pj + 3*sg + 1*rg + i] + x[pj + 3*sg + 2*rg + i] +
                           x[pj + 4*sg + 0*rg + i] + x[pj + 4*sg + 1*rg + i] + x[pj + 4*sg + 2*rg + i] +
                           x[pj + 5*sg + 0*rg + i] + x[pj + 5*sg + 1*rg + i] + x[pj + 5*sg + 2*rg + i])
          }
          
          for (i in 22:24) {
            PS[i] <- sum(x[pj + 3*sg + r*rg + i] + x[pj + 4*sg + r*rg + i] + x[pj + 5*sg + r*rg + i])
          }
        } else {
          # Rcpp::Rcout << "Place where a >= 12." 
          for (i in 1:24) {
            PS[i] <- sum(x[pj + 3*sg + 0*rg + i] + x[pj + 3*sg + 1*rg + i] + x[pj + 3*sg + 2*rg + i] +
                           x[pj + 4*sg + 0*rg + i] + x[pj + 4*sg + 1*rg + i] + x[pj + 4*sg + 2*rg + i] +
                           x[pj + 5*sg + 0*rg + i] + x[pj + 5*sg + 1*rg + i] + x[pj + 5*sg + 2*rg + i])
          }
        }
##        
        cpmu <- cpo <- cMmu <- cMo <- cLo <- cMa <- cpmu_dose <- cpo_dose <- cMmu_dose <- cMo_dose <- cLo_dose <- cMa_dose <- 0
        
        p <- a * ag + s * sg + r * rg
        q <- a * ag + kpc + r * rg
        o <- ifelse(s > 3, a * ag + (s - 3) * sg + r * rg, 0)
        
        if (s < 3) {
          cpmu <- 0
          cpmu_dose <- 0
          cpo <- 0
          cpo_dose <- 0
          cMmu <- 0
          cMmu_dose <- 0
          cMo <- 0
          cMo_dose <- 0
          cLo <- 0
          cLo_dose <- 0
          cMa <- 0
          cMa_dose <- 0
          SN <- 1
          lossMS1 <- lossMS0 <- lossP <- 0
        } else {
          lossP <- x[p + 21] * (1.0 / 60.0)
          lossMS0 <- x[p + 22] * om_mab
          lossMS1 <- x[p + 23] * om_mab
          
          cpmu <- 0
          cpmu_dose <- 0
          cpo <- 0
          cpo_dose <- 0
          cMmu <- 0
          cMmu_dose <- 0
          cMo <- 0
          cMo_dose <- 0
          cLo <- 0
          cLo_dose <- 0
          cMa <- 0
          cMa_dose <- 0
          
          if (a == 0) {
            if (r == 0) {
              cpo <- cal_pal(t1, 0)
              cMo <- cal_mAB_VHR(t1, 0)
              cpo_dose <- cal_pal_dose(t1, 0)
              cMo_dose <- cal_mAB_VHR_dose(t1, 0)
            } else if (r == 1) {
              cMo <- cal_mAB_HR(t1, 0)
              cMo_dose <- cal_mAB_HR_dose(t1, 0)
            } else {
              cMo <- cal_mAB_LR(t1, 0)
              cMo_dose <- cal_mAB_LR_dose(t1, 0)
            }
          } else if (a > 0) {
            if (r == 0) {
              cpo <- cal_pal(t1, a)
              cpo_dose <- cal_pal_dose(t1, a)
              cMo <- cal_mAB_VHR(t1, a)
              cMo_dose <- cal_mAB_VHR_dose(t1, a)
              S0 <- x[a * ag + (s - 3) * sg + r * rg + 1]
              S1 <- x[a * ag + (s - 3) * sg + r * rg + 6]
              S2 <- x[a * ag + (s - 3) * sg + r * rg + 11]
              S3 <- x[a * ag + (s - 3) * sg + r * rg + 16]
              SN <- 1
            } else if (r == 1) {
              cMo <- cal_mAB_HR(t1, a)
              cLo <- cal_LAV_HR(t1, a)
              cMo_dose <- cal_mAB_HR_dose(t1, a)
              cLo_dose <- cal_LAV_HR_dose(t1, a)
              S0 <- x[a * ag + (s - 3) * sg + r * rg + 1]
              S1 <- x[a * ag + (s - 3) * sg + r * rg + 6]
              S2 <- x[a * ag + (s - 3) * sg + r * rg + 11]
              S3 <- x[a * ag + (s - 3) * sg + r * rg + 16]
              SN <- 1
            } else {
              cMo <- cal_mAB_LR(t1, a)
              cMo_dose <- cal_mAB_LR_dose(t1, a)
              if (s == 4) {
                cMa <- vac_cal(t1, a)
                cMa_dose <- vac_cal_dose(t1, a)
                S0 <- x[a * ag + (s - 3) * sg + r * rg + 1]
                S1 <- x[a * ag + (s - 3) * sg + r * rg + 6]
                S2 <- x[a * ag + (s - 3) * sg + r * rg + 11]
                S3 <- x[a * ag + (s - 3) * sg + r * rg + 16]
              } else {
                cMa <- 0
                cMa_dose <- 0
                SN <- 1
                cLo <- cal_LAV_LR(t1, a)
                cLo_dose <- cal_LAV_LR_dose(t1, a)
                S0 <- x[a * ag + (s - 3) * sg + r * rg + 1]
                S1 <- x[a * ag + (s - 3) * sg + r * rg + 6]
                S2 <- x[a * ag + (s - 3) * sg + r * rg + 11]
                S3 <- x[a * ag + (s - 3) * sg + r * rg + 16]
              }
            }
          }
        else
        {
            
        }
        }
        for (i in 1:21) {
          x_tot <- x_tot + x[o + i]
          PST <- PST + PS[i]
        }
        
        if (PST < 1) {
          PST <- 1
        }
        
        x_tot_1 <- x_tot_2 <- 0
        
        for (i in 1:2) {
          x_tot_1 <- x_tot_1 + x[o + i]
        }
        
        for (i in 2:5) {
          x_tot_2 <- x_tot_2 + x[o + i]
        }
        
        # Keeping track of doses
        protectpal <- protectpal + cpo_dose * x_tot
        protectmabs <- protectmabs + cMo_dose * x_tot
        protectLAV <- protectLAV + cLo_dose * x_tot
        protectmat <- protectmat + cMa_dose * x_tot
        
        # Differential Equations
        dxdt[p + 0] <- (1.0 - p_vul) * mu * rp + mu_mat * rp - x[p + 0] * xi * xi_b - x[p + 0] * ej1 + PS[0] * ej * rp * u - x[o + 0] * cMo - x[o + 0] * cpo
        dxdt[p + 1] <- p_vul * mu * rp + x[p + 0] * xi * xi_b + lossP + lossMS0 - x[p + 1] * In * beta - x[p + 1] * ej1 + PS[1] * ej * rp * u - x[o + 1] * cMo - x[o + 1] * cpo - cLo * S0 - cMa * S0
        dxdt[p + 2] <- x[p + 1] * In * beta - x[p + 2] * si - x[p + 2] * ej1 + PS[2] * ej * rp * u - x[o + 2] * cMo - x[o + 2] * cpo
        dxdt[p + 3] <- x[p + 2] * si * pA[a] - x[p + 3] * ga0 * rho - x[p + 3] * ej1 + PS[3] * ej * rp * u - x[o + 3] * cMo - x[o + 3] * cpo
        dxdt[p + 4] <- x[p + 2] * si * (1.0 - pA[a]) - x[p + 4] * ga0 - x[p + 4] * ej1 + PS[4] * ej * rp * u - x[o + 4] * cMo - x[o + 4] * cpo
        dxdt[p + 5] <- x[p + 4] * ga0 + x[p + 3] * ga0 * rho - x[p + 5] * om - x[p + 5] * ej1 + PS[5] * ej * rp * u - x[o + 5] * cMo - x[o + 5] * cpo + cLo * S0 + cMa * S0
        dxdt[p + 6] <- x[p + 5] * om - d1 * x[p + 6] * In * beta + lossMS1 - x[p + 6] * ej1 + PS[6] * ej * rp * u - x[o + 6] * cMo - x[o + 6] * cpo - cLo * S1 - cMa * S1
        dxdt[p + 7] <- d1 * x[p + 6] * In * beta - x[p + 7] * si - x[p + 7] * ej1 + PS[7] * ej * rp * u - x[o + 7] * 0 - x[o + 7] * cpo
        dxdt[p + 8] <- x[p + 7] * si * pA[a] - x[p + 8] * ga1 * rho - x[p + 8] * ej1 + PS[8] * ej * rp * u - x[o + 8] * 0 - x[o + 8] * cpo
        dxdt[p + 9] <- x[p + 7] * si * (1.0 - pA[a]) - x[p + 9] * ga1 - x[p + 9] * ej1 + PS[9] * ej * rp * u - x[o + 9] * 0 - x[o + 9] * cpo
        dxdt[p + 10] <- x[p + 9] * ga1 + x[p + 8] * ga1 * rho - x[p + 10] * om - x[p + 10] * ej1 + PS[10] * ej * rp * u - x[o + 10] * 0 - x[o + 10] * cpo + cLo * S1 + cMa * S1
        dxdt[p + 11] <- x[p + 10] * om - d2 * x[p + 11] * In * beta - x[p + 11] * ej1 + PS[11] * ej * rp * u - x[o + 11] * 0 - x[o + 11] * cpo - cLo * S2 - cMa * S2
        dxdt[p + 12] <- d2 * x[p + 11] * In * beta - x[p + 12] * si - x[p + 12] * ej1 + PS[12] * ej * rp * u - x[o + 12] * 0 - x[o + 12] * cpo
        dxdt[p + 13] <- x[p + 12] * si * pA[a] - x[p + 13] * ga2 * rho - x[p + 13] * ej1 + PS[13] * ej * rp * u - x[o + 13] * 0 - x[o + 13] * cpo
        dxdt[p + 14] <- x[p + 12] * si * (1.0 - pA[a]) - x[p + 14] * ga2 - x[p + 14] * ej1 + PS[14] * ej * rp * u - x[o + 14] * 0 - x[o + 14] * cpo
        dxdt[p + 15] <- x[p + 14] * ga2 + x[p + 13] * ga2 * rho - x[p + 15] * om - x[p + 15] * ej1 + PS[15] * ej * rp * u - x[o + 15] * 0 - x[o + 15] * cpo + cLo * S2 + cMa * S2
        dxdt[p + 16] <- x[p + 15] * om + x[p + 20] * om - d3 * x[p + 16] * In * beta - x[p + 16] * ej1 + PS[16] * ej * rp * u - x[o + 16] * 0 - x[o + 16] * cpo - cLo * S3 - cMa * S3
        dxdt[p + 17] <- d3 * x[p + 16] * In * beta - x[p + 17] * si - x[p + 17] * ej1 + PS[17] * ej * rp * u - x[o + 17] * 0 - x[o + 17] * cpo
        dxdt[p + 18] <- x[p + 17] * si * pA[a] - x[p + 18] * ga3 * rho - x[p + 18] * ej1 + PS[18] * ej * rp * u - x[o + 18] * 0 - x[o + 18] * cpo
        dxdt[p + 19] <- x[p + 17] * si * (1.0 - pA[a]) - x[p + 19] * ga3 - x[p + 19] * ej1 + PS[19] * ej * rp * u - x[o + 19] * 0 - x[o + 19] * cpo
        dxdt[p + 20] <- x[p + 19] * ga3 + x[p + 18] * ga3 * rho - x[p + 20] * om - x[p + 20] * ej1 + PS[20] * ej * rp * u - x[o + 20] * 0 - x[o + 20] * cpo + cLo * S3 + cMa * S3
        
        # Vaccination Groups
        dxdt[p + 21] <- x_tot * cpo - lossP - x[p + 21] * ej1 + PS[21] * ej * u * rp  # pal
        dxdt[p + 22] <- x_tot_1 * cMo - lossMS0 - x[p + 22] * ej1 + PS[22] * ej * u * rp  # mabs
        dxdt[p + 23] <- x_tot_2 * cMo - lossMS1 - x[p + 23] * ej1 + PS[23] * ej * u * rp  # mabs
  }
}
# Vaccine groups
dxdt[a * ag + 6 * sg + 0] <- si * (x[cj + 3 * sg + 0 * rg + 2] + x[cj + 4 * sg + 0 * rg + 2] + x[cj + 5 * sg + 0 * rg + 2] + x[cj + 3 * sg + 1 * rg + 2] + x[cj + 4 * sg + 1 * rg + 2] + x[cj + 5 * sg + 1 * rg + 2] + x[cj + 3 * sg + 2 * rg + 2] + x[cj + 4 * sg + 2 * rg + 2] + x[cj + 5 * sg + 2 * rg + 2])
dxdt[a * ag + 6 * sg + 1] <- protectpal
dxdt[a * ag + 6 * sg + 2] <- protectmabs
dxdt[a * ag + 6 * sg + 3] <- protectLAV
dxdt[a * ag + 6 * sg + 4] <- protectmat
dxdt[a * ag + 6 * sg + 5] <- 0
dxdt[a * ag + 6 * sg + 6] <- 0
dxdt[a * ag + 6 * sg + 7] <- 0

# Monitoring Parent
dxdt[a * ag + 6 * sg + 8] <- si * (x[cj + 3 * sg + 0 * rg + 2] + x[cj + 3 * sg + 0 * rg + 7] + x[cj + 3 * sg + 0 * rg + 12] + x[cj + 3 * sg + 0 * rg + 17])  # VHR
dxdt[a * ag + 6 * sg + 9] <- si * (x[cj + 3 * sg + 1 * rg + 2] + x[cj + 3 * sg + 1 * rg + 7] + x[cj + 3 * sg + 1 * rg + 12] + x[cj + 3 * sg + 1 * rg + 17])  # HR
dxdt[a * ag + 6 * sg + 10] <- si * (x[cj + 3 * sg + 2 * rg + 2] + x[cj + 3 * sg + 2 * rg + 7] + x[cj + 3 * sg + 2 * rg + 12] + x[cj + 3 * sg + 2 * rg + 17])  # LR

# Monitoring Cocoon
dxdt[a * ag + 6 * sg + 11] <- si * (x[cj + 4 * sg + 0 * rg + 2] + x[cj + 4 * sg + 0 * rg + 7] + x[cj + 4 * sg + 0 * rg + 12] + x[cj + 4 * sg + 0 * rg + 17])  # VHR
dxdt[a * ag + 6 * sg + 12] <- si * (x[cj + 4 * sg + 1 * rg + 2] + x[cj + 4 * sg + 1 * rg + 7] + x[cj + 4 * sg + 1 * rg + 12] + x[cj + 4 * sg + 1 * rg + 17])  # HR
dxdt[a * ag + 6 * sg + 13] <- si * (x[cj + 4 * sg + 2 * rg + 2] + x[cj + 4 * sg + 2 * rg + 7] + x[cj + 4 * sg + 2 * rg + 12] + x[cj + 4 * sg + 2 * rg + 17])  # LR

# Monitoring Neither
dxdt[a * ag + 6 * sg + 14] <- si * (x[cj + 5 * sg + 0 * rg + 2] + x[cj + 5 * sg + 0 * rg + 7] + x[cj + 5 * sg + 0 * rg + 12] + x[cj + 5 * sg + 0 * rg + 17])  # VHR
dxdt[a * ag + 6 * sg + 15] <- si * (x[cj + 5 * sg + 1 * rg + 2] + x[cj + 5 * sg + 1 * rg + 7] + x[cj + 5 * sg + 1 * rg + 12] + x[cj + 5 * sg + 1 * rg + 17])  # HR
dxdt[a * ag + 6 * sg + 16] <- si * (x[cj + 5 * sg + 2 * rg + 2] + x[cj + 5 * sg + 2 * rg + 7] + x[cj + 5 * sg + 2 * rg + 12] + x[cj + 5 * sg + 2 * rg + 17])  # LR

x[a * ag + 6 * sg + 17] <- N_tot_n[a]
x[a * ag + 6 * sg + 18] <- N_tot_c[a]
x[a * ag + 6 * sg + 19] <- N_tot_p[a]

x[a * ag + 6 * sg + 20] <- N_tot_n_v[a]
x[a * ag + 6 * sg + 21] <- N_tot_c_v[a]
x[a * ag + 6 * sg + 22] <- N_tot_p_v[a]
   } 
  }
 }
}