rm(list = ls())
library(deSolve)
############################ parameters #####################
params <- list(b1           = 1.998,                            # relative amplitude
               phi          = 0.614,                           # seasonal offset/peak transmission seasonal
               psi          = 0.236,                            # width of seasonal peak/standard deviation
               sigma        = 1/4.98,                          # susceptibility to primary infection
               d0           = 1,                                  # relative susceptibility to primary infection [d0=1 not included in Hodgson et al, added for uniformity]
               d1           = 0.89,                           # relative susceptibility to secondary infection, relative to  primary infection
               d2           = 0.81,                           # relative susceptibility to tertiary infection, relative to secondary infection
               d3           = 0.33,                           # relative susceptibility to subsequent infections after third infection, relative to tertiary infection
               alpha        = 0.634,                           # reduction in infectiousness of asymptomatic infections
               propR        = 0.5,                                # proportion of neonates born with protection
               mu           = 1863,                                # daily birth rate # 1/1863
               xi           = 1/133.5,                          # rate of loss maternal immunity
               probA        = c(0.255, 0.635,0.753),                                # probability infection is asymptomatic
               g0           = 6.062,                             # rate of loss of primary infectiousness
               g1           = 0.878,                           # proportional decrease between secondary and primary infection
               g2           = 0.812,                           # proportional decrease between tertiary and secondary infection
               gamma3       = 1,                                  # duration of infectiousness for all infections after tertiary infection
               omega        = 1/358.9,                          # rate of loss of disease-induced/post-infection immunity
               age_r_1        = 1/(365*4),
               age_r_2        = 1/(365*59),
               age_r_3        = 1/(365*10),           # aging rate 
               waning_rate1 = 0,                        # vaccine-induced waning rate age group 1
               waning_rate2 = 0,                        # vaccine-induced waning rate age group 2
               waning_rate3 = 0,          # vaccine-induced waning rate age group 3
               uptake1 = 0,                               # uptake rate of intervention age group 1
               uptake2 = 0,                               # uptake rate of intervention age group 2
               uptake3 = 0,                               # uptake rate of intervention age group 3
               vac_eff1 = 0,                              # vaccine efficacy age group 1
               vac_eff2 = 0,                              # vaccine efficacy age group 2
               vac_eff3 = 0      # vaccine efficacy age group 3
               
)
############################ function to calculate gamma #####################

calculate_gamma(params = params)
############################ age adjustment initial states##################### 
# Define the age ranges
ranges <- list(c(1, 16), c(17, 23), c(24, 25))

# function to sum across ranges
sum_ranges <- function(data, ranges) {
  sums <- lapply(ranges, function(range) {
    age_group_sum <- rep(0, length(data[[1]]))
    
    for (i in range[1]:range[2]) {
      age_group_sum <- age_group_sum + unlist(data[i])
    }
    
    return(age_group_sum)
  })
  return(sums)
}

# merged initial states 
state_age_merged <- sum_ranges(state_initial, ranges)

# check for consistency
sum(initial_states)              # matrix original 25 age groups = 55744600
sum(unlist(state_initial))   # list original 25 age groups = 55744600
sum(unlist(state_age_merged))    # list merged to 3 groups     = 55744600
sum(round(unlist(state_age_merged), digits = 0))  # rounded sum = 55744594 
# manually input states_initial_age from 'state_age_merged' above
state_initial_age <- c(M=c(244302,0,0),    S0=c(318282,34100,0),  E0=c(200915,21526,0),   A0=c(23465,9013,0),  I0=c(137742,8259,0), R0=c(360445,38618,0), V0=c(0,0,0),
                       S1=c(302992,126300,0),  E1=c(151606,63196,0), A1=c(22426,30270,0),  I1=c(116178,27507,0), R1=c(343130,143031,0), V1=c(0,0,0),
                       S2=c(217519,261480,0),   E2=c(70910,85241,0), A2=c(13360,50703,0),  I2=c(66456,45244,0),   R2=c(246333,296118,0), V2=c(0,0,0),
                       S3=c(240637,16972110,4979010), E3=c(22164,1563208,458589),   A3=c(4218,1243715,380711), I3=c(20730 ,515845,135481),   R3=c(272514,19220418,5638583), V3=c(0,0,0),
                       Z=c(0,0,0))

# state_initial_age <- c(M=c(244302,0,0),    S0=c(349583,2799,0),  E0=c(220674,1767,0),   A0=c(31715,763,0),  I0=c(145346,655,0), R0=c(395893,3170,0), V0=c(0,0,0),
#                        S1=c(411519,17773,0),  E1=c(205909,8893,0), A1=c(48260,4436,0),  I1=c(139991,3694,0), R1=c(466034,20127,0), V1=c(0,0,0),
#                        S2=c(419062,59936,0),   E2=c(136612,19539,0), A2=c(51843,12220,0),  I2=c(101928,9773,0),   R2=c(474576,67876,0), V2=c(0,0,0),
#                        S3=c(1266790,15945956,4979010), E3=c(116677,1468694,458589),   A3=c(59576,1188357,380711), I3=c(71757,464819,135481),   R3=c(1434603,18058329,5638583), V3=c(0,0,0),
#                        Z=c(0,0,0))

############################ function to update_ages #####################
update_ages <- function(pop,age_r){
  age_pop = pop*age_r
  dP = - age_pop +  c(0,age_pop[-length(age_pop)])
  P = pop + dP
  return(list(P=P,dP=dP))
}
############################ rsv model age updated function #####################
rsv_model_age <- function(t, state_initial_age, params){
  #  with(as.list(c(state_initial,params)),{
  M  = state_initial_age[grepl('M',names(state_initial_age))]
  
  S0 = state_initial_age[grepl('S0',names(state_initial_age))]
  E0 = state_initial_age[grepl('E0',names(state_initial_age))]
  A0 = state_initial_age[grepl('A0',names(state_initial_age))]
  I0 = state_initial_age[grepl('I0',names(state_initial_age))]
  R0 = state_initial_age[grepl('R0',names(state_initial_age))]
  V0 = state_initial_age[grepl('V0',names(state_initial_age))]
  
  S1 = state_initial_age[grepl('S1',names(state_initial_age))]
  E1 = state_initial_age[grepl('E1',names(state_initial_age))]
  A1 = state_initial_age[grepl('A1',names(state_initial_age))]
  I1 = state_initial_age[grepl('I1',names(state_initial_age))]
  R1 = state_initial_age[grepl('R1',names(state_initial_age))]
  V1 = state_initial_age[grepl('V1',names(state_initial_age))]
  
  S2 = state_initial_age[grepl('S2',names(state_initial_age))]
  E2 = state_initial_age[grepl('E2',names(state_initial_age))]
  A2 = state_initial_age[grepl('A2',names(state_initial_age))]
  I2 = state_initial_age[grepl('I2',names(state_initial_age))]
  R2 = state_initial_age[grepl('R2',names(state_initial_age))]
  V2 = state_initial_age[grepl('V2',names(state_initial_age))]
  
  S3 = state_initial_age[grepl('S3',names(state_initial_age))]
  E3 = state_initial_age[grepl('E3',names(state_initial_age))]
  A3 = state_initial_age[grepl('A3',names(state_initial_age))]
  I3 = state_initial_age[grepl('I3',names(state_initial_age))]
  R3 = state_initial_age[grepl('R3',names(state_initial_age))]
  V3 = state_initial_age[grepl('V3',names(state_initial_age))]
  
  Z  = state_initial_age[grepl('Z',names(state_initial_age))]
  
  # parameters 
  propR  = params[['propR']]
  mu     = params[['mu']]
  xi     = params[['xi']]
  probA  = params[['probA']]
  alpha  = params[['alpha']]
  sigma  = params[['sigma']]
  omega  = params[['omega']]
  gamma0 = params[['gamma0']] 
  gamma1 = params[['gamma1']] 
  gamma2 = params[['gamma2']] 
  gamma3 = params[['gamma3']]
  b1     = params[['b1']]
  phi    = params[['phi']]
  psi    = params[['psi']]
  d0     = params[['d0']]
  d1     = params[['d1']]
  d2     = params[['d2']]
  d3     = params[['d3']]
  qp     = params[['qp']]
  qc     = params[['qc']]
  age_r  = c(params[['age_r_1']], params[['age_r_2']], params[['age_r_3']])
  N      = sum(unlist(state_initial_age))
  uptake = c(params[['uptake1']], params[['uptake2']], params[['uptake3']]) 
  vac_eff = c(params[['vac_eff1']], params[['vac_eff2']], params[['vac_eff3']]) 
  waning_rate = c(params[['waning_rate1']], params[['waning_rate2']], params[['waning_rate3']]) 
  
  # age update initial states 
  ageM = update_ages(M,age_r)
  dM   = ageM$dP
  M    = ageM$P
  
  ageS0 = update_ages(S0,age_r)
  dS0   = ageS0$dP
  S0    = ageS0$P
  ageE0 = update_ages(E0,age_r)
  dE0   = ageE0$dP
  E0    = ageE0$P
  ageA0 = update_ages(A0,age_r)
  dA0   = ageA0$dP
  A0    = ageA0$P
  ageI0 = update_ages(I0,age_r)
  dI0   = ageI0$dP
  I0    = ageI0$P
  ageR0 = update_ages(R0,age_r)
  dR0   = ageR0$dP
  R0    = ageR0$P
  ageV0 = update_ages(V0,age_r)
  dV0   = ageV0$dP
  V0    = ageV0$P
  
  ageS1 = update_ages(S1,age_r)
  dS1   = ageS1$dP
  S1    = ageS1$P
  ageE1 = update_ages(E1,age_r)
  dE1   = ageE1$dP
  E1    = ageE1$P
  ageA1 = update_ages(A1,age_r)
  dA1   = ageA1$dP
  A1    = ageA1$P
  ageI1 = update_ages(I1,age_r)
  dI1   = ageI1$dP
  I1    = ageI1$P
  ageR1 = update_ages(R1,age_r)
  dR1   = ageR1$dP
  R1    = ageR1$P
  ageV1 = update_ages(V1,age_r)
  dV1   = ageV1$dP
  V1    = ageV1$P
  
  ageS2 = update_ages(S2,age_r)
  dS2   = ageS2$dP
  S2    = ageS2$P
  ageE2 = update_ages(E2,age_r)
  dE2   = ageE2$dP
  E2    = ageE2$P
  ageA2 = update_ages(A2,age_r)
  dA2   = ageA2$dP
  A2    = ageA2$P
  ageI2 = update_ages(I2,age_r)
  dI2   = ageI2$dP
  I2    = ageI2$P
  ageR2 = update_ages(R2,age_r)
  dR2   = ageR2$dP
  R2    = ageR2$P
  ageV2 = update_ages(V2,age_r)
  dV2   = ageV2$dP
  V2    = ageV2$P
  
  ageS3 = update_ages(S3,age_r)
  dS3   = ageS3$dP
  S3    = ageS3$P
  ageE3 = update_ages(E3,age_r)
  dE3   = ageE3$dP
  E3    = ageE3$P
  ageA3 = update_ages(A3,age_r)
  dA3   = ageA3$dP
  A3    = ageA3$P
  ageI3 = update_ages(I3,age_r)
  dI3   = ageI3$dP
  I3    = ageI3$P
  ageR3 = update_ages(R3,age_r)
  dR3   = ageR3$dP
  R3    = ageR3$P
  ageV3 = update_ages(V3,age_r)
  dV3   = ageV3$dP
  V3    = ageV3$P
  
  # vaccination
  # get day of the year
  t1 = t %% 365
  
  year_length = 365
  vac_day = 300 
  start_year = 6  
  end_year = 10   
  
  # start and end day for vaccination
  start_vaccination_day = vac_day + (start_year * year_length)
  end_vaccination_day = vac_day + (end_year - 1) * year_length
  
  if (t >= start_vaccination_day && t <= end_vaccination_day && abs(t1 - vac_day) <= 1) {

    # dV0_vac = ((S0 + E0 + A0) * uptake* vac_eff)
    # # V0 = V0 + dV0_vac 
    # dV0 = dV0 + dV0_vac
    dS0_vac = (S0 * uptake* vac_eff)
    # S0 = S0 - dS0_vac 
    dS0 = dS0 - dS0_vac
    dE0_vac = (E0 * uptake* vac_eff)
    # E0 = E0 - dE0_vac
    dE0 = dE0 - dE0_vac
    dA0_vac = (A0 * uptake* vac_eff)
    # A0 = A0 - dA0_vac
    dA0 = dA0 - dA0_vac
    dV0_vac = dV0 + dS0_vac + dE0_vac + dA0_vac
    dV0 = dV0 + dV0_vac
    
    # dV1_vac = ((S1 + E1 + A1) * uptake* vac_eff)
    # # V1 = V1 + dV1_vac 
    # dV1 = dV1 + dV1_vac
    dS1_vac = (S1 * uptake* vac_eff)
    # S1 = S1 - dS1_vac 
    dS1 = dS1 - dS1_vac
    dE1_vac = (E1 * uptake* vac_eff)
    # E1 = E1 - dE1_vac
    dE1 = dE1 - dE1_vac
    dA1_vac = (A1 * uptake* vac_eff)
    # A1 = A1 - dA1_vac
    dA1 = dA1 - dA1_vac
    dV1_vac = dV1 + dS1_vac + dE1_vac + dA1_vac
    dV1 = dV1 + dV1_vac

    # dV2_vac = ((S2 + E2 + A2) * uptake* vac_eff)
    # # V2 = V2 + dV2_vac 
    # dV2 = dV2 + dV2_vac
    dS2_vac = (S2 * uptake* vac_eff)
    # S2 = S2 - dS2_vac 
    dS2 = dS2 - dS2_vac
    dE2_vac = (E2 * uptake* vac_eff)
    # E2 = E2 - dE0_vac
    dE2 = dE2 - dE2_vac
    dA2_vac = (A2 * uptake* vac_eff)
    # A2 = A2 - dA2_vac
    dA2 = dA2 - dA2_vac
    dV2_vac = dV2 + dS2_vac + dE2_vac + dA2_vac
    dV2 = dV2 + dV2_vac

    # # dV3_vac = ((S3 + E3 + A3) * uptake* vac_eff)
    # # # V3 = V3 + dV3_vac
    # # dV3 = dV3 + dV3_vac
    # dS3_vac = (S3 * uptake* vac_eff)
    # # S3 = S3 + dS3_vac
    # dS3 = dS3 - dS3_vac
    # dE3_vac = (E3 * uptake* vac_eff)
    # # E3 = E3 + dE3_vac
    # dE3 = dE3 - dE3_vac
    # dA3_vac = (A3 * uptake* vac_eff)
    # # A3 = A3 + dA3_vac
    # dA3 = dA3 - dA3_vac
    # dV3_vac = dV3 + dS3_vac + dE3_vac + dA3_vac
    # dV3 = dV3 + dV3_vac
  }
  
  # beta
  beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))

  # lambda
  lambda0 = beta_value * d0 *(sum(I0)+(sum(A0)*alpha)) / sum(S0,E0,I0,A0,R0)
  lambda1 = beta_value * d1* (sum(I1)+(sum(A1)*alpha)) / sum(S1,E1,I1,A1,R1)
  lambda2 = beta_value * (d1*d2)*(sum(I2)+(sum(A2)*alpha)) / sum(S2,E2,I2,A2,R2)
  lambda3 = beta_value * (d1*d2*d3)*(sum(I3)+(sum(A3)*alpha)) / sum(S3,E3,I3,A3,R3)
  
  
  # age adjusted health - ODEs
  
  threshold = 1
  
  dM  = dM + propR*mu*c(1,0,0) - M*xi
  
  dS0 = dS0 + (1-propR)*mu*c(1,0,0) + M*xi - S0*lambda0  + V0*waning_rate
  dE0 = dE0 + S0*lambda0 - E0*sigma 
  dA0 = dA0 + E0*probA*sigma - A0*gamma0 
  dI0 = dI0 + E0*(1-probA)*sigma - I0*gamma0 
  dR0 = dR0 + I0*gamma0 + A0*gamma0 - R0*omega 
  dV0 = dV0 - V0*waning_rate 
  
  dS1 = dS1 + R0*omega - S1*lambda1 + V1*waning_rate
  dE1 = dE1 + S1*lambda1 - E1*sigma 
  dA1 = dA1 + E1*probA*sigma - A1*gamma1 
  dI1 = dI1 + E1*(1-probA)*sigma - I1*gamma1
  dR1 = dR1 + I1*gamma1 + A1*gamma1 - R1*omega
  dV1 = dV1 - V1*waning_rate 
  
  dS2 = dS2 + R1*omega - S2*lambda2  + V2*waning_rate
  dE2 = dE2 + S2*lambda2 - E2*sigma 
  dA2 = dA2 + E2*probA*sigma - A2*gamma2
  dI2 = dI2 + E2*(1-probA)*sigma - I2*gamma2
  dR2 = dR2 + I2*gamma2 + A2*gamma2 - R2*omega
  dV2 = dV2 - V2*waning_rate 
  
  dS3 = dS3 + R2*omega + R3*omega - S3*lambda3 + V3*waning_rate
  dE3 = dE3 + dE3 + S3*lambda3 - E3*sigma 
  dA3 = dA3 + E3*probA*sigma - A3*gamma3 
  dI3 = dI3 + E3*(1-probA)*sigma - I3*gamma3
  dR3 = dR3 + I3*gamma3 + A3*gamma3 - R3*omega
  dV3 = dV3 - V3*waning_rate 
  
  dZ = sigma * (E0+E1+E2+E3)
  
  return(list(c(dM,dS0,dE0,dA0,dI0,dR0,dV0,dS1,dE1,dA1,dI1,dR1,dV1,dS2,dE2,dA2,dI2,dR2,dV2,dS3,dE3,dA3,dI3,dR3,dV3,dZ)))
  #})
}
