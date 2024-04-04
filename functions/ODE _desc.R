### source functions 
source("/Users/vincentarumadri/Desktop/Epi/Modelling/master_thesis_ua/functions/RunInterventions.R")

#### ODE start
ODE_desc <- function(vac_calendar, vac_dose, vac_info, cov_c) {
  
# Perform calculations and assign additional properties
xi <- 1.0 / parameterValues[["xi"]]
si <- 1.0 / parameterValues[["si"]]
ga0 <- 1.0 / parameterValues[["ga0"]]
ga1 <- 1.0 / (parameterValues[["ga0"]] * parameterValues[["g1"]])
ga2 <- 1.0 / (parameterValues[["ga0"]] * parameterValues[["g1"]] * parameterValues[["g2"]])
ga3 <- ga2
om <- 1.0 / parameterValues[["om"]]
  
rho <- 1.0
alpha_i <- parameterValues[["alpha_i"]]
d1 <- parameterValues[["d1"]]
d2 <- parameterValues[["d1"]] * parameterValues[["d2"]]
d3 <- parameterValues[["d1"]] * parameterValues[["d2"]] * parameterValues[["d3"]]
a1 <- 1.0
a2 <- 1.0
a3 <- 1.0
  
phi <- parameterValues[["phi"]]
qp <- parameterValues[["qp"]]
qc <- parameterValues[["qc"]]
b1 <- parameterValues[["b1"]]
psi <- parameterValues[["psi"]]
  
A <- A
eta <- eta
  
u18p <- u_p[1]
u19p <- u_p[2]
u20p <- u_p[3]
  
populationPerAgeGroup <- populationPerAgeGroup
pVHR <- pVHR
pHR <- pHR
pLR <- pLR
p_mat <- p_mat

dailyBirthRate <- dailyBirthRate
pA <- pA
  
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
  
PS <- rg
  
cal_mAB_VHR_t <- vac_calendar[["mAB_VHR"]]
cal_mAB_LR_t <- vac_calendar[["mAB_LR"]]
cal_LAV_LR_t <- vac_calendar[["LAV_LR"]]
cal_mat_LR_t <- vac_calendar[["mat_LR"]]
  
cal_mAB_VHR_step_out_t <- vac_calendar[["mAB_VHR_custom_out"]]
cal_mAB_LR_step_out_t <- vac_calendar[["mAB_LR_custom_out"]]
  
cal_mAB_VHR_dose <- vac_dose[["mAB_VHR"]]
cal_mAB_LR_dose <- vac_dose[["mAB_LR"]]
cal_LAV_LR_dose <- vac_dose[["LAV_LR"]]
vac_cal_dose <- vac_dose[["mat_LR"]]
  
cal_mAB_VHR_step_out_dose <- vac_dose[["mAB_VHR_custom_out"]]
cal_mAB_LR_step_out_dose <- vac_dose[["mAB_LR_custom_out"]]
  
cal_mAB_HR_step_out_dose <- matrix(0, nrow = 365, ncol = A)
cal_mAB_HR_dose <- matrix(0, nrow = 365, ncol = A)
cal_LAV_HR_dose <- matrix(0, nrow = 365, ncol = A)
  
direct <- vac_info[["direct"]]
  
# Sorting out efficacy values
eff_wane_lav <- vac_info[["lav_mass"]]
eff_wane_mat <- vac_info[["mat_mass"]]
eff_wane_mab <- vac_info[["mab_mass"]]
eff_wane_vhr <- vac_info[["mab_vhr"]]
eff_wane_step <- 0
  
xi_boost <- 1
  
rg <- rg
sg <- sg
ag <- ag
}  
######## operator    
model_rsv <- function(t, state, parms) {
    with(as.list(c(state, parms)), {
    # CALENDARS
    vac_cal_vhr <- matrix(0, 365, A)  
    vac_cal <- matrix(0, 365, A)
    
    if (t < 2 * 365) {
      cal_pal <- matrix(0, 365, A)
      vac_cal <- matrix(0, 365, A)
      cal_mAB_VHR <- matrix(0, 365, A)
      cal_mAB_HR <- matrix(0, 365, A)
      cal_mAB_LR <- matrix(0, 365, A)
      cal_LAV_HR <- matrix(0, 365, A)
      cal_LAV_LR <- matrix(0, 365, A)
      
      # THE STEP FUNCTION IMMUNITY 
      cal_mAB_VHR_step_out <- matrix(0, nrow = 365, ncol = A)
      cal_mAB_HR_step_out <- matrix(0, nrow = 365, ncol = A)
      cal_mAB_LR_step_out <- matrix(0, nrow = 365, ncol = A)
      
    } else {
      cal_mAB_VHR <- cal_mAB_VHR_t
      cal_mAB_LR <- cal_mAB_LR_t
      cal_LAV_LR <- cal_LAV_LR_t
      vac_cal <- cal_mat_LR_t
      
      # THE STEP FUNCTION IMMUNITY
      cal_mAB_VHR_step_out <- cal_mAB_VHR_step_out_t
      cal_mAB_LR_step_out <- cal_mAB_LR_step_out_t
      
    }

##### populations
tot_coc <- 0
tot_coc_v <- 0
tot_inf <- 0
tot_inf_v <- 0

for (a in 1:A) {
  
  if (a < 12) {  # Maybe be 13 to adjust for R's 1-based indexing ## need to check 
      N_tot_n[a] <- populationPerAgeGroup[a] * (1 - cov_c)
      N_tot_c[a] <- populationPerAgeGroup[a] * cov_c
      N_tot_p[a] <- 0
      
      N_tot_n_v[a] <- populationPerAgeGroup[a] * (1 - cov_c)
      N_tot_c_v[a] <- populationPerAgeGroup[a] * cov_c
      N_tot_p_v[a] <- 0
    } else {
      N_tot_n[a] <- populationPerAgeGroup[a] * (1 - p_mat[a])
      N_tot_c[a] <- populationPerAgeGroup[a] * p_mat[a] * cov_c
      N_tot_p[a] <- populationPerAgeGroup[a] * p_mat[a] * (1 - cov_c)
      
      N_tot_n_v[a] <- populationPerAgeGroup[a] * (1 - p_mat[a])
      N_tot_c_v[a] <- populationPerAgeGroup[a] * p_mat[a] * cov_c
      N_tot_p_v[a] <- populationPerAgeGroup[a] * p_mat[a] * (1 - cov_c)
  }
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

# Initialization
protectpal <- 0
protectmabs <- 0
protectLAV <- 0
protectmat <- 0
step_in <- step_out <- step_out_dose <- step_in_S0 <- step_in_S1 <- 0
loss_mat <- loss_lav_S0 <- loss_lav_S1 <- loss_lav_S2 <- loss_lav_S3 <- 0

for (a in 1:A) {
  r_prop <- 0
  protectpal <- protectmabs <- protectLAV <- protectmat <- 0 <- step_out_dose <- step_in <- step_out <- step_in_S0 <- step_in_S1 <- 0
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
  lossMS2 <- 0
  
  # Birth rate into each social group
  if (a == 1) {
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
        u <- 1
        In <- I_temp_c
        mu <- muBc
        mu_mat <- 0 
        p_vul <- p_vulc 
        xi_b <-  xi_bc
      } else if (s == 3) {
        kp <- 2 * sg
        u <- 1
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
        u <- 1
        In <- I_temp_c_v
        mu <- muBcv * (1 - vac_cal(t1, 0))
        mu_mat <- muBcv*vac_cal(t1,0)
        p_vul <- p_vulcv 
        xi_b <-  xi_bcv
      } else {
        kp <- 5 * sg
        u <- 1 
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
      # "Risk groups parameters" 
      if ((a == 1) & (s == 5)) {   
        kp <- 4 * sg 
        if (r == 1) {
          mu <- muBcv 
          mu_mat <- 0
        } else {
          mu <- muBcv * (1 - vac_cal[t1, 1]) 
          mu_mat <- muBcv * vac_cal[t1, 1]
        }
      } else if ((a == 1) & (s == 6)) {
        kp <- 5 * sg 
        if (r == 1) {
          mu <- muBnv 
          mu_mat <- 0
        } else {
          mu <- muBnv 
          mu_mat <- 0
        }
      }
      
      PST <- 0
      x_tot <- x_tot_1 <- x_tot_2 <- x_tot_3 <- 0
      step_out <- step_in <- step_in_S0 <- step_in_S1 <- 0
      
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
        #"Place where r == 3." 
        if (a %in% c(18, 19, 20)) {
          up <- switch(u18p, u19p, u20p)
        }
      }
      
      if (s < 3) {
        # "Place where s < 3." 
        if (a < 12) {
          for (i in 1:21) {
            PS[i] <- x[pj + s*sg + 0*rg + i] + x[pj + s*sg + 1*rg + i] + x[pj + s*sg + 2*rg + i] 
          }
          
          for (i in 22:rg) {
            PS[i] <- x[pj + s*sg + 0*rg + i] + x[pj + s*sg + 1*rg + i] + x[pj + s*sg + 2*rg + i]
          }
        } else {
          for (i in 1:rg) {
            PS[i] <- (x[pj + 0*sg + 0*rg + i] + x[pj + 0*sg + 1*rg + i] + x[pj + 0*sg + 2*rg + i] +
                           x[pj + 1*sg + 0*rg + i] + x[pj + 1*sg + 1*rg + i] + x[pj + 1*sg + 2*rg + i] +
                           x[pj + 2*sg + 0*rg + i] + x[pj + 2*sg + 1*rg + i] + x[pj + 2*sg + 2*rg + i])
          }
        }
      } else {
        # "Place where s >= 3." 
        if (a < 12) {
          # "Place where a < 12." 
          for (i in 1:21) {
            PS[i] <- x[pj + s*sg + 0*rg + i] + x[pj + s*sg + 1*rg + i] + x[pj + s*sg + 2*rg + i] 
          
          }
          for (i in 22:rg) {
            PS[i] <- sum(x[pj + s*sg + 0*rg + i] + x[pj + s*sg + 1*rg + i] + x[pj + s*sg + 2*rg + i])
          }
        } else {
          # "Place where a >= 12." 
          for (i in 1:rg) {
            PS[i] <- sum(x[pj + 3*sg + 0*rg + i] + x[pj + 3*sg + 1*rg + i] + x[pj + 3*sg + 2*rg + i] +
                           x[pj + 4*sg + 0*rg + i] + x[pj + 4*sg + 1*rg + i] + x[pj + 4*sg + 2*rg + i] +
                           x[pj + 5*sg + 0*rg + i] + x[pj + 5*sg + 1*rg + i] + x[pj + 5*sg + 2*rg + i])
          }
        }
          }
## p defined        
cpmu <- cpo <- cMmu <- cMo <- cLo <- cMa <- cpmu_dose <- cpo_dose <- cMmu_dose <- cMo_dose <- cLo_dose <- cMa_dose <- 0
p <- a * ag + s * sg + r * rg
q <- a * ag + kpc + r * rg
o <- ifelse(s > 3, a * ag + (s - 3) * sg + r * rg, 0)

if (s > 3) 
   o = a * ag + (s-3) * sg + r * rg
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
          lossMS2 <- lossMS1 <- lossMS0 <- lossP <- 0
          step_in <- step_out <- step_in_S0 <- step_in_S1 <- 0
          loss_mat <- loss_lav_S0 <-  loss_lav_S1 <- loss_lav_S2 <- loss_lav_S3 <- lossMS0 <- lossMS1 <- lossMS2 <- 0
        } else {
          if (r == 1) {
            loss_mat <- loss_lav_S0 <- loss_lav_S1 <- loss_lav_S2 <- loss_lav_S3 <- 0
            lossMS0 <- x[p + 26] * eff_wane_vhr
            lossMS1 <- x[p + 29] * eff_wane_vhr
            lossMS2 <- x[p + 32] * eff_wane_vhr
          } else if (r == 3) {
            loss_mat <- x[p + 23] * eff_wane_mat
            lossMS0 <- x[p + 26] * eff_wane_mab
            lossMS1 <- x[p + 29] * eff_wane_mab
            lossMS2 <- x[p + 32] * eff_wane_mab
            loss_lav_S0 <- x[p + 35] * eff_wane_lav
            loss_lav_S1 <- x[p + 38] * eff_wane_lav
            loss_lav_S2 <- x[p + 41] * eff_wane_lav
            loss_lav_S3 <- x[p + 44] * eff_wane_lav
          } else if (r == 2) {
            loss_mat <- loss_lav_S0 <- loss_lav_S1 <- loss_lav_S2 <- loss_lav_S3 <- lossMS0 <- lossMS1 <- lossMS2 <- 0
          }
          
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
step_in <- step_out <- step_in_S0 <- step_in_S1 <- 0

          if (a == 1) {
            M <- x[a*ag + (s-3)*sg + r*rg + 0]
            S0 <- x[a*ag + (s-3)*sg + r*rg + 1]
            S1 <- x[a*ag + (s-3)*sg + r*rg + 6]
            S2 <- x[a*ag + (s-3)*sg + r*rg + 11]
            S3 <- x[a*ag + (s-3)*sg + r*rg + 16]
            
            if (r == 1) {
              if(s == 5) {
                step_in_S0 <- 0 
                step_in_S1 <- 0 
                step_out_dose <- step_out_dose + cal_mAB_VHR_step_out_dose[t1, 1] # Adding 1 to R index
              }
              step_out <- cal_mAB_VHR_step_out[t1, 1] 
              
              cMo <- cal_mAB_VHR(t1, 1)
              cMo_dose <- cal_mAB_VHR_dose(t1, 1)
            } else if (r == 2) {
              if(s == 5) {
                step_in_S0 <- 0 
                step_in_S1 <- 0 
                step_out_dose <- step_out_dose + cal_mAB_HR_step_out_dose[t1, 1] 
              }
              step_out <- cal_mAB_HR_step_out[t1, 1] 
              cMo <- cal_mAB_HR(t1, 1)
              cMo_dose <- cal_mAB_HR_dose(t1, 1)
            } else { 
              if(s == 5) {
                step_in_S0 <- collect_protect(t, 4) + collect_protect(t, 0) + collect_protect(t, 2)
                step_in_S1 <- collect_protect(t, 5) + collect_protect(t, 1) + collect_protect(t, 3)
                step_out_dose <- step_out_dose + cal_mAB_LR_step_out_dose[t1, 1] 
              }
              step_out <- cal_mAB_LR_step_out[t1, 1] 
              cMo <- cal_mAB_LR(t1, 0)
              cMo_dose <- cal_mAB_LR_dose(t1, 0)
            }
          } else if (a > 1) {
            M <- x[a*ag + (s-3)*sg + r*rg + 0]
            if (r == 1) {
              if(s == 5) {
                step_in_S0 = 0 # collect_protect(t, 6 * (a) + 0)
                step_in_S1 = 0 # collect_protect(t, 6 * (a) + 1)
                step_out_dose = step_out_dose + cal_mAB_VHR_step_out_dose[t1, a + 1] 
              }
              step_out = cal_mAB_VHR_step_out[t1, a + 1] 
              
              cpo <- cal_pal(t1, a)
              cpo_dose <- cal_pal_dose(t1, a)
              cMo <- cal_mAB_VHR(t1, a)
              cMo_dose <- cal_mAB_VHR_dose(t1, a)
              S0 <- x[a * ag + (s - 3) * sg + r * rg + 1]
              S1 <- x[a * ag + (s - 3) * sg + r * rg + 6]
              S2 <- x[a * ag + (s - 3) * sg + r * rg + 11]
              S3 <- x[a * ag + (s - 3) * sg + r * rg + 16]
              SN <- 1
            } else if (r == 2) {
              if(s == 5) {
                step_in_S0 = 0 # collect_protect(t, 6 * (a) + 2)
                step_in_S1 = 0 # collect_protect(t, 6 * (a) + 3)
                step_out_dose = step_out_dose + cal_mAB_HR_step_out_dose[t1, a + 1] 
              }
              step_out = cal_mAB_HR_step_out[t1, a + 1] 
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
              if(s == 5) {
                step_in_S0 = collect_protect(t, 6 * a + 4 + 1) + collect_protect(t, 6 * a + 0 + 1) + collect_protect(t, 6 * a + 2 + 1) 
                step_in_S1 = collect_protect(t, 6 * a + 5 + 1) + collect_protect(t, 6 * a + 1 + 1) + collect_protect(t, 6 * a + 3 + 1) 
                step_out_dose = step_out_dose + cal_mAB_LR_step_out_dose[t1, a + 1] 
              }
              step_out = cal_mAB_LR_step_out[t1, a + 1] 
              cMo <- cal_mAB_LR(t1, a)
              cMo_dose <- cal_mAB_LR_dose(t1, a)
              
              if (s == 4) {
                cLo <- vac_cal(t1, a)
                cLo_dose <- vac_cal_dose(t1, a)
                S0 <- x[a * ag + (s - 3) * sg + r * rg + 1]
                S1 <- x[a * ag + (s - 3) * sg + r * rg + 6]
                S2 <- x[a * ag + (s - 3) * sg + r * rg + 11]
                S3 <- x[a * ag + (s - 3) * sg + r * rg + 16]
              } else {
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
        
# Keeping track of doses
protectpal <- protectpal + cpo_dose * x_tot
protectmabs <- protectmabs + cMo_dose * x_tot
protectLAV <- protectLAV + cLo_dose * x_tot
protectmat <- protectmat + cMa_dose * x_tot
            
  if (s == 5) {
              if (r == 1) {
                collect_protect[t + eff_wane_step, 6 * (a + 5) + 0] <- x_tot_1 * step_out
                collect_protect[t + eff_wane_step, 6 * (a + 5) + 1] <- x_tot_2 * step_out
              }
              
            if (r == 2) {
              collect_protect(t + eff_wane_step, 6 * (a + 5) + 2) <- x_tot_1 * step_out
              collect_protect(t + eff_wane_step, 6 * (a + 5) + 3) <- x_tot_2 * step_out
            } else if (r == 3) {
              collect_protect(t + eff_wane_step, 6 * (a + 5) + 4) <- x_tot_1 * step_out
              collect_protect(t + eff_wane_step, 6 * (a + 5) + 5) <- x_tot_2 * step_out
            }
          }
            mab_trans <- numeric(12)
            lav_trans <- numeric(4)
            
            for (i in 1:12) {
              state <- x[p + i - 1] 
              prot <- x[o + i - 1] * cMo 
              if ((state - prot) < 0) {
                mab_trans[i] <- x[p + i - 1] * cMo 
              } else {
                mab_trans[i] <- prot
              }
            }
            
            for (i in 1:2) {
              x_tot_1 <- x_tot_1 + mab_trans[i]
            }
            
            for (i in 3:7) {
              x_tot_2 <- x_tot_2 + mab_trans[i]
            }
            
            for (i in 8:12) {
              x_tot_3 <- x_tot_3 + mab_trans[i]
            }
            
            pos <- c(2, 7, 12, 17)
            
            for (i in 1:4) {
              state <- x[p + pos[i]]
              prot <- x[o + pos[i]] * cLo
              if ((state - prot) < 0) {
                lav_trans[i] <- x[p + pos[i]] * cLo
              } else {
                lav_trans[i] <- prot
              }
            }
            
            # Differential Equations
            dxdt[p+0] = (1.0-p_vul)*mu*rp - x[p+0]*xi*xi_b - (x[p+0])*ej1 + PS[0]*ej*rp*u - mab_trans[0]
            
            dxdt[p+1] = p_vul*mu*rp  + x[p+0]*xi*xi_b + lossMS0 - x[p+1]*In*beta - (x[p+1])*ej1 + PS[1]*ej*rp*u - mab_trans[1] - x[o+1]*step_out - lav_trans[0] + step_in_S0 + x[p+23]*(eff_wane_mat)
            dxdt[p+2] = x[p+1]*In*beta                 - x[p+2]*si          - (x[p+2])*ej1 + PS[2]*ej*rp*u -mab_trans[2] - x[p+2]*step_out
            dxdt[p+3] = x[p+2]*si*pA[a]                - x[p+3]*ga0*rho     -  (x[p+3])*ej1 + PS[3]*ej*rp*u - mab_trans[3] - x[p+3]*step_out
            dxdt[p+4] = x[p+2]*si*(1.0-pA[a])          - x[p+4]*ga0         - (x[p+4])*ej1 + PS[4]*ej*rp*u - mab_trans[4] - x[p+4]*step_out
            dxdt[p+5] = x[p+4]*ga0 + x[p+3]*ga0*rho    - x[p+5]*om          - (x[p+5])*ej1 + PS[5]*ej*rp*u - mab_trans[5] - x[p+5]*step_out
            
            dxdt[p+6] = x[p+5]*om                      - d1*x[p+6]*In*beta + lossMS1  - (x[p+6])*ej1 + PS[6]*ej*rp*u - mab_trans[6] + loss_lav_S0 - lav_trans[1] + step_in_S1
            dxdt[p+7] = d1*x[p+6]*In*beta              - x[p+7]*si          - (x[p+7])*ej1 + PS[7]*ej*rp*u - mab_trans[7]
            dxdt[p+8] = x[p+7]*si*pA[a]                - x[p+8]*ga1*rho     - (x[p+8])*ej1 + PS[8]*ej*rp*u - mab_trans[8]
            dxdt[p+9] = x[p+7]*si*(1.0-pA[a])          - x[p+9]*ga1         - (x[p+9])*ej1 + PS[9]*ej*rp*u - mab_trans[9]
            dxdt[p+10] = x[p+9]*ga1 + x[p+8]*ga1*rho   - x[p+10]*om         - (x[p+10])*ej1 + PS[10]*ej*rp*u - mab_trans[10]
            
            dxdt[p+11] = x[p+10]*om                    - d2*x[p+11]*In*beta + lossMS2 - (x[p+11])*ej1 + PS[11]*ej*rp*u - mab_trans[11] + loss_lav_S1 - lav_trans[2]
            dxdt[p+12] = d2*x[p+11]*In*beta            - x[p+12]*si         - (x[p+12])*ej1 + PS[12]*ej*rp*u
            dxdt[p+13] = x[p+12]*si*pA[a]              - x[p+13]*ga2*rho    - (x[p+13])*ej1 + PS[13]*ej*rp*u
            dxdt[p+14] = x[p+12]*si*(1.0-pA[a])        - x[p+14]*ga2        - (x[p+14])*ej1 + PS[14]*ej*rp*u
            dxdt[p+15] = x[p+14]*ga2 + x[p+13]*ga2*rho - x[p+15]*om         - (x[p+15])*ej1 + PS[15]*ej*rp*u
            
            dxdt[p+16] = x[p+15]*om + x[p+20]*om       - d3*x[p+16]*In*beta - (x[p+16])*ej1 + PS[16]*ej*rp*u  + loss_lav_S2 + loss_lav_S3 - lav_trans[3]
            dxdt[p+17] = d3*x[p+16]*In*beta            - x[p+17]*si         - (x[p+17])*ej1 + PS[17]*ej*rp*u
            dxdt[p+18] = x[p+17]*si*pA[a]              - x[p+18]*ga3*rho    - (x[p+18])*ej1 + PS[18]*ej*rp*u
            dxdt[p+19] = x[p+17]*si*(1.0-pA[a])        - x[p+19]*ga3        - (x[p+19])*ej1 + PS[19]*ej*rp*u
            dxdt[p+20] = x[p+19]*ga3 + x[p+18]*ga3*rho - x[p+20]*om         - (x[p+20])*ej1 + PS[20]*ej*rp*u
            
            # Maternal vaccination
            dxdt[p+21] = mu_mat*rp - x[p+21]*eff_wane_mat - x[p+21]*ej1 + PS[21]*ej*u*rp
            dxdt[p+22] = x[p+21]*eff_wane_mat - x[p+22]*eff_wane_mat - x[p+22]*ej1 + PS[22]*ej*u*rp
            dxdt[p+23] = x[p+22]*eff_wane_mat - x[p+23]*eff_wane_mat - x[p+23]*ej1 + PS[23]*ej*u*rp
            
            # Monoclonal protection
            dxdt[p+24] = x_tot_1 - x[p+24]*eff_wane_mab  - x[p+24]*ej1 + PS[24]*ej*u*rp
            dxdt[p+25] = x[p+24]*eff_wane_mab - x[p+25]*eff_wane_mab  - x[p+25]*ej1 + PS[25]*ej*u*rp
            dxdt[p+26] = x[p+25]*eff_wane_mab - lossMS0 - x[p+26]*ej1 + PS[26]*ej*u*rp
            
            dxdt[p+27] = x_tot_2 - x[p+27]*eff_wane_mab - x[p+27]*ej1 + PS[27]*ej*u*rp
            dxdt[p+28] = x[p+27]*eff_wane_mab - x[p+28]*eff_wane_mab - x[p+28]*ej1 + PS[28]*ej*u*rp
            dxdt[p+29] = x[p+28]*eff_wane_mab - lossMS1 - x[p+29]*ej1 + PS[29]*ej*u*rp
            
            dxdt[p+30] = x_tot_3 - x[p+30]*eff_wane_mab - x[p+30]*ej1 + PS[30]*ej*u*rp
            dxdt[p+31] = x[p+30]*eff_wane_mab - x[p+31]*eff_wane_mab - x[p+31]*ej1 + PS[31]*ej*u*rp
            dxdt[p+32] = x[p+31]*eff_wane_mab - lossMS2 - x[p+32]*ej1 + PS[32]*ej*u*rp
            
            # LAV protection
            dxdt[p+33] = lav_trans[0] - x[p+33]*eff_wane_lav - x[p+33]*ej1 + PS[33]*ej*u*rp
            dxdt[p+34] = x[p+33]*eff_wane_lav - x[p+34]*eff_wane_lav - x[p+34]*ej1 + PS[34]*ej*u*rp
            dxdt[p+35] = x[p+34]*eff_wane_lav - loss_lav_S0 - x[p+35]*ej1 + PS[35]*ej*u*rp
            
            dxdt[p+36] = lav_trans[1] - x[p+36]*eff_wane_lav - x[p+36]*ej1 + PS[36]*ej*u*rp
            dxdt[p+37] = x[p+36]*eff_wane_lav - x[p+37]*eff_wane_lav - x[p+37]*ej1 + PS[37]*ej*u*rp
            dxdt[p+38] = x[p+37]*eff_wane_lav - loss_lav_S1 - x[p+38]*ej1 + PS[38]*ej*u*rp
            
            dxdt[p+39] = lav_trans[2] - x[p+39]*eff_wane_lav - x[p+39]*ej1 + PS[39]*ej*u*rp
            dxdt[p+40] = x[p+39]*eff_wane_lav - x[p+40]*eff_wane_lav - x[p+40]*ej1 + PS[40]*ej*u*rp
            dxdt[p+41] = x[p+40]*eff_wane_lav - loss_lav_S2 - x[p+41]*ej1 + PS[41]*ej*u*rp
            
            dxdt[p+42] = lav_trans[3] - x[p+42]*eff_wane_lav - x[p+42]*ej1 + PS[42]*ej*u*rp
            dxdt[p+43] = x[p+42]*eff_wane_lav - x[p+43]*eff_wane_lav - x[p+43]*ej1 + PS[43]*ej*u*rp
            dxdt[p+44] = x[p+43]*eff_wane_lav - loss_lav_S3 - x[p+44]*ej1 + PS[44]*ej*u*rp
        } 
      }
            # Vaccine groups
            dxdt[a*ag + 6*sg + 0] = si*(x[cj+3*sg+0*rg+2] + x[cj+4*sg+0*rg+2] + x[cj+5*sg+0*rg+2] + x[cj+3*sg+1*rg+2] + x[cj+4*sg+1*rg+2] + x[cj+5*sg+1*rg+2] + x[cj+3*sg+2*rg+2] + x[cj+4*sg+2*rg+2] + x[cj+5*sg+2*rg+2])
            dxdt[a*ag + 6*sg + 1] = protectpal
            dxdt[a*ag + 6*sg + 2] = protectmabs + step_out_dose
            dxdt[a*ag + 6*sg + 3] = protectLAV
            
            dxdt[a*ag + 6*sg + 4] = protectmat
            dxdt[a*ag + 6*sg + 5] = 0
            dxdt[a*ag + 6*sg + 6] = 0
            dxdt[a*ag + 6*sg + 7] = 0
            
            # Monitoring Parent
            dxdt[a*ag + 6*sg + 8] =  si*(x[cj+3*sg+0*rg+2]+x[cj+3*sg+0*rg+7]+x[cj+3*sg+0*rg+12]+x[cj+3*sg+0*rg+17]) # VHR
            dxdt[a*ag + 6*sg + 9] =  si*(x[cj+3*sg+1*rg+2]+x[cj+3*sg+1*rg+7]+x[cj+3*sg+1*rg+12]+x[cj+3*sg+1*rg+17]) # HR
            dxdt[a*ag + 6*sg + 10] = si*(x[cj+3*sg+2*rg+2]+x[cj+3*sg+2*rg+7]+x[cj+3*sg+2*rg+12]+x[cj+3*sg+2*rg+17]) # LR
      
            # Monitoring Cocoon
            dxdt[a*ag + 6*sg + 11] = si*(x[cj+4*sg+0*rg+2]+x[cj+4*sg+0*rg+7]+x[cj+4*sg+0*rg+12]+x[cj+4*sg+0*rg+17]) # VHR
            dxdt[a*ag + 6*sg + 12] = si*(x[cj+4*sg+1*rg+2]+x[cj+4*sg+1*rg+7]+x[cj+4*sg+1*rg+12]+x[cj+4*sg+1*rg+17]) # HR
            dxdt[a*ag + 6*sg + 13] = si*(x[cj+4*sg+2*rg+2]+x[cj+4*sg+2*rg+7]+x[cj+4*sg+2*rg+12]+x[cj+4*sg+2*rg+17]) # LR
      
            # Monitoring Neither
            dxdt[a*ag + 6*sg + 14] = si*(x[cj+5*sg+0*rg+2]+x[cj+5*sg+0*rg+7]+x[cj+5*sg+0*rg+12]+x[cj+5*sg+0*rg+17]) # VHR
            dxdt[a*ag + 6*sg + 15] = si*(x[cj+5*sg+1*rg+2]+x[cj+5*sg+1*rg+7]+x[cj+5*sg+1*rg+12]+x[cj+5*sg+1*rg+17]) # HR
            dxdt[a*ag + 6*sg + 16] = si*(x[cj+5*sg+2*rg+2]+x[cj+5*sg+2*rg+7]+x[cj+5*sg+2*rg+12]+x[cj+5*sg+2*rg+17]) # LR
      
            x[a*ag + 6*sg + 17] = N_tot_n[a]
            x[a*ag + 6*sg + 18] = N_tot_c[a]
            x[a*ag + 6*sg + 19] = N_tot_p[a]
      
            x[a*ag + 6*sg + 20] = N_tot_n_v[a]
            x[a*ag + 6*sg + 21] = N_tot_c_v[a]
            x[a*ag + 6*sg + 22] = N_tot_p_v[a]
            
          }
 # Return the list of derivatives
  list(dxdt)
    })
}
###
check_stability <- function(x0) {
  neg_indices <- which(x0 < 0)
  
  if (length(neg_indices) > 0) {
    for (i in neg_indices) {
      cat(sprintf("i: %d. x0: %f\n", i, x0[i]))
      x0[i] = 0
    }
  }
  
  return(x0)
}
