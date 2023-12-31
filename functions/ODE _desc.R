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
#### ODE function ####
ODE_desc <- function(finELL_t, vac_calendar_t, vac_dose_t, vac_info_t, cov_c_t) {
  
  source("functions/Run_Interventions.R")
  
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

