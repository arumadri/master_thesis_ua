rm(list = ls())
############################ parameters #####################
params_check <- list(b1     = 1.998,                            # relative amplitude
               phi          = 0.614,                            # seasonal offset/peak transmission seasonal
               psi          = 0.236,                            # width of seasonal peak/standard deviation
               sigma        = 1/4.98,                           # susceptibility to primary infection
               d0           = 1,                                # relative susceptibility to primary infection [d0=1 not included in Hodgson et al, added for uniformity]
               d1           = 0.89,                             # relative susceptibility to secondary infection, relative to  primary infection
               d2           = 0.81,                             # relative susceptibility to tertiary infection, relative to secondary infection
               d3           = 0.33,                             # relative susceptibility to subsequent infections after third infection, relative to tertiary infection
               alpha        = 0.634,                            # reduction in infectiousness of asymptomatic infections
               propR        = 0.5,                              # proportion of neonates born with protection
               mu           = 1863,                             # daily birth rate 
               xi           = 1/133.5,                          # rate of loss maternal immunity
               probA        = c(0.127, 0.635,0.753),            # probability infection is asymptomatic
               g0           = 6.16,                             # rate of loss of primary infectiousness
               g1           = 0.87,                             # decrease in secondary infection duration relative to primary
               g2           = 0.79,                             # decrease in subsequent infection duration relative to secondary
               gamma3       = 1,                                # duration of infectiousness for all infections after tertiary infection
               omega        = 1/358.9,                          # rate of loss of disease-induced/post-infection immunity
               age_r_1      = 1/(365*5),                        # aging rate, age group 1
               age_r_2      = 1/(365*60),                       # aging rate, age group 2
               age_r_3      = 1/(365*25),                       # aging rate, age group 3
               waning_rate1 = 0,                                # vaccine-induced waning rate, age group 1
               waning_rate2 = 0,                                # vaccine-induced waning rate, age group 2
               waning_rate3 = 0,                                # vaccine-induced waning rate, age group 3
               uptake1      = 0.8,                              # uptake rate of intervention, age group 1
               uptake2      = 0.8,                              # uptake rate of intervention, age group 2
               uptake3      = 0.8,                              # uptake rate of intervention, age group 3
               vac_eff1     = 0.6,                              # vaccine efficacy, age group 1
               vac_eff2     = 0.6,                              # vaccine efficacy, age group 2
               vac_eff3     = 0.6                               # vaccine efficacy, age group 3
               
)
calculate_gamma <- function(params) {
  
  # parameters from list
  g0 = params_check[['g0']]
  g1 = params_check[['g1']]
  g2 = params_check[['g2']]
  
  # new parameters
  params_check[['gamma0']] <<- 1/g0
  params_check[['gamma1']] <<- 1/(g0*g1)
  params_check[['gamma2']] <<- 1/(g0*g1*g2)
  
}
calculate_gamma(params = params_check)

############################ age adjustment initial states##################### 
state_initial_check <- c(M=c(100,50,20), S0=c(100,50,20), E0=c(100,50,20), A0=c(100,50,20), I0=c(100,50,20),  R0=c(100,50,20), V0 = c(0,0,0),
                       S1=c(100,50,20), E1=c(100,50,20), A1=c(100,50,20), I1=c(100,50,20),  R1=c(100,50,20), V1 = c(0,0,0),
                       S2=c(100,50,20), E2=c(100,50,20), A2=c(100,50,20), I2=c(100,50,20),  R2=c(100,50,20), V2 = c(0,0,0),
                       S3=c(100,50,20), E3=c(100,50,20), A3=c(100,50,20), I3=c(100,50,20),  R3=c(100,50,20), V3 = c(0,0,0),
                       Z=c(0,0,0))
############################ function to update_ages #####################
update_ages <- function(pop,age_r){
  age_pop = pop*age_r
  dP = - age_pop +  c(0,age_pop[-length(age_pop)])
  P = pop + dP
  return(list(P=P,dP=dP))
}
############################ rsv model age updated function #####################
rsv_model_check <- function(t, state_initial_check, params){
  #  with(as.list(c(state_initial,params)),{
  M  = state_initial_check[grepl('M',names(state_initial_check))]
  
  S0 = state_initial_check[grepl('S0',names(state_initial_check))]
  E0 = state_initial_check[grepl('E0',names(state_initial_check))]
  A0 = state_initial_check[grepl('A0',names(state_initial_check))]
  I0 = state_initial_check[grepl('I0',names(state_initial_check))]
  R0 = state_initial_check[grepl('R0',names(state_initial_check))]
  V0 = state_initial_check[grepl('V0',names(state_initial_check))]
  
  S1 = state_initial_check[grepl('S1',names(state_initial_check))]
  E1 = state_initial_check[grepl('E1',names(state_initial_check))]
  A1 = state_initial_check[grepl('A1',names(state_initial_check))]
  I1 = state_initial_check[grepl('I1',names(state_initial_check))]
  R1 = state_initial_check[grepl('R1',names(state_initial_check))]
  V1 = state_initial_check[grepl('V1',names(state_initial_check))]
  
  S2 = state_initial_check[grepl('S2',names(state_initial_check))]
  E2 = state_initial_check[grepl('E2',names(state_initial_check))]
  A2 = state_initial_check[grepl('A2',names(state_initial_check))]
  I2 = state_initial_check[grepl('I2',names(state_initial_check))]
  R2 = state_initial_check[grepl('R2',names(state_initial_check))]
  V2 = state_initial_check[grepl('V2',names(state_initial_check))]
  
  S3 = state_initial_check[grepl('S3',names(state_initial_check))]
  E3 = state_initial_check[grepl('E3',names(state_initial_check))]
  A3 = state_initial_check[grepl('A3',names(state_initial_check))]
  I3 = state_initial_check[grepl('I3',names(state_initial_check))]
  R3 = state_initial_check[grepl('R3',names(state_initial_check))]
  V3 = state_initial_check[grepl('V3',names(state_initial_check))]
  
  Z  = state_initial_check[grepl('Z',names(state_initial_check))]
  
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
  age_r  = c(params[['age_r_1']], params[['age_r_2']], params[['age_r_3']])
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
  
  ageZ = update_ages(Z,age_r)
  dZ   = ageZ$dP
  Z    = ageZ$P
  
  # vaccination
  # get day of the year
  t1 = t %% 365

  year_length = 365
  vac_day = 274  # 1st October, start of RSV season in UK, assumption is that eligible individuals are vaccinated by this date
  start_year = 0
  end_year = 1  # intervention for one season

  # start and end day for vaccination
  start_vaccination_day = vac_day + (start_year * year_length)
  end_vaccination_day = (end_year * year_length)

  if (t >= start_vaccination_day && t <= end_vaccination_day && abs(t1 - vac_day) <= 1) {

    dS0_vac = (S0 * uptake* vac_eff)
    dS0     = dS0 - dS0_vac
    dV0     = dV0 + dS0_vac

    dS1_vac = (S1 * uptake* vac_eff)
    dS1     = dS1 - dS1_vac
    dV1     = dV1 + dS1_vac

    dS2_vac = (S2 * uptake* vac_eff)
    dS2     = dS2 - dS2_vac
    dV2     = dV2 + dS2_vac

    # dS3_vac = (S3 * uptake* vac_eff)  ## exclude vaccination in exposure level 3
    # dS3     = dS3 - dS3_vac
    # dV3     = dV3 + dS3_vac
  }
  
  # beta
  beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))
  
  # lambda (foi)
  lambda0 = beta_value * d0 *(sum(I0)+(sum(A0)*alpha)) / sum(S0,E0,I0,A0,R0)
  lambda1 = beta_value * d1* (sum(I1)+(sum(A1)*alpha)) / sum(S1,E1,I1,A1,R1)
  lambda2 = beta_value * (d1*d2)*(sum(I2)+(sum(A2)*alpha)) / sum(S2,E2,I2,A2,R2)
  lambda3 = beta_value * (d1*d2*d3)*(sum(I3)+(sum(A3)*alpha)) / sum(S3,E3,I3,A3,R3)
  
  
  # age adjusted health ODEs
  
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
  dE3 = dE3 + S3*lambda3 - E3*sigma 
  dA3 = dA3 + E3*probA*sigma - A3*gamma3 
  dI3 = dI3 + E3*(1-probA)*sigma - I3*gamma3
  dR3 = dR3 + I3*gamma3 + A3*gamma3 - R3*omega
  dV3 = dV3 - V3*waning_rate 
  
  dZ = sigma * (E0+E1+E2+E3)
  
  return(list(c(dM,dS0,dE0,dA0,dI0,dR0,dV0,dS1,dE1,dA1,dI1,dR1,dV1,dS2,dE2,dA2,dI2,dR2,dV2,dS3,dE3,dA3,dI3,dR3,dV3,dZ)))
  #})
}

## 1. checking vaccination
# health parameters = 0, aging_rate = 0, only change intervention (t=274)
params_check[1:18] <- 0
params_check$age_r_1 <- 0
params_check$age_r_2 <- 0
params_check$age_r_3 <- 0
params_check$d0 <- 0
params_check$gamma0 <- 0
params_check$gamma1 <- 0
params_check$gamma2 <- 0
params_check$gamma3 <- 0

# check vaccination implementation
matrix(unlist(rsv_model_check(274,state_initial_check,params_check)),ncol=3,byrow = T) 

# above = movement out of susceptible compartment (2,8,14) into vaccinated compartments (7,13,19), no vaccination in exposure level 3, hence no movement 

## 2. checking health 
# To check health, comment out the vaccination implementation code in model, line 198:223, set aging = 0
params_check <- params
params_check$age_r_1 <- 0
params_check$age_r_2 <- 0
params_check$age_r_3 <- 0

# check health implementation
matrix(unlist(rsv_model_check(365,state_initial_check,params_check)),ncol=3,byrow = T) 

# above = only health changes, no aging, no vaccination hence no movement in vaccination compartments (7,13,19,25)
