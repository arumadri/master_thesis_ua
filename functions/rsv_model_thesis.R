rm(list = ls())
############################ parameters #####################
params <- list( 
               b1           = 1.997732,                            # mean seasonal transmission rate
               phi          = 0.615044,                           # seasonal offset/peak transmission seasonal
               psi          = 0.2401791,                            # width of seasonal peak/standard deviation
               sigma        = 1/4.863965,                          # susceptibility to primary infection
               d0           = 1,                                  # relative susceptibility to primary infection [d0=1 not included in Hodgson et al, added for uniformity]
               d1           = 0.8908731,                           # relative susceptibility to secondary infection, relative to  primary infection
               d2           = 0.7912774,                           # relative susceptibility to tertiary infection, relative to secondary infection
               d3           = 0.3428649,                           # relative susceptibility to subsequent infections after third infection, relative to tertiary infection
               alpha        = 0.6085629,                           # reduction in infectiousness of asymptomatic infections
               propR        = 0.3,                                # proportion of neonates born with protection
               mu           = 1/1863,                                # daily birth rate # 1/1863
               xi           = 1/131.2745,                          # rate of loss maternal immunity
               probA        = 0.5,                                # probability infection is asymptomatic
               g0           = 6.062045,                             # rate of loss of primary infectiousness
               g1           = 0.8776286,                           # proportional decrease between secondary and primary infection
               g2           = 0.8122186,                           # proportional decrease between tertiary and secondary infection
               gamma3       = 1,                                  # duration of infectiousness for all infections after tertiary infection
               omega        = 1/356.8321,                          # rate of loss of disease-induced/post-infection immunity
               age_r        = 1/365,                              # aging rate
               waning_rate = 1/180,                              # intervention waning rate - age group 1, generic will be updated in scenario analysis
               # waning_rate2 = 1/180,                              # intervention waning rate - age group 2, generic will be updated in scenario analysis
               # waning_rate3 = 1/180,                              # intervention waning rate - age group 3, generic will be updated in scenario analysis
               uptake      = 0,                                # uptake rate of intervention - age group 1, generic will be updated in scenario analysis
               # uptake2      = 0,                                # uptake rate of intervention - age group 2, generic will be updated in scenario analysis
               # uptake3      = 0,                                # uptake rate of intervention - age group 3, generic will be updated in scenario analysis
               eff         = 0                               # intervention efficacy - age group 1, generic will be updated in scenario analysis
               # eff2         = 0,                               # intervention efficacy - age group 2, generic will be updated in scenario analysis
               # eff3         = 0                                # intervention efficacy - age group 3, generic will be updated in scenario analysis
) 
############################ function to calculate gamma #####################
calculate_gamma <- function(params) {
  
  # parameters from list
  g0 = params[['g0']]
  g1 = params[['g1']]
  g2 = params[['g2']]
  
  # store new parameters
  params[['gamma0']] <<- 1/g0
  params[['gamma1']] <<- 1/(g0*g1)
  params[['gamma2']] <<- 1/(g0*g1*g2)
  
}
calculate_gamma(params)
# ############################ age adjustment initial states #####################
# state_initial_age <- c(M=c(244302,0,0), S0=c(349583,34100,0), E0=c(220674,21526,0), A0=c(31715,9013,0), I0=c(145346,8259,0),  R0=c(395893,38618,0), V0 = c(0,0,0),
#                                         S1=c(411519,126300,0), E1=c(205909,63196,0), A1=c(48260,30270,0), I1=c(139991,27507,0),  R1=c(466034,143031,0), V1 = c(0,0,0),
#                                         S2=c(419062,261480,0), E2=c(136612,85241,0), A2=c(51843,50703,0), I2=c(101928,45244,0),  R2=c(474576,296118,0), V2 = c(0,0,0),
#                                         S3=c(1266790,16972110,7896503), E3=c(116677,1563208,727303), A3=c(59576,1243715,603793), I3=c(71757,515845,214867),  R3=c(1434603,19220418,8942559), V3 = c(0,0,0),
#                        Z=c(0,0,0))
# ############################ function to update_ages #####################
update_ages <- function(pop,age_r){
  age_pop = pop*age_r
  dP = - age_pop +  c(0,age_pop[-length(age_pop)])
  P = pop + dP
  return(list(P=P,dP=dP))
}
############################ rsv model age updated function #####################
rsv_model_age <- function(t, state_initial_age, params){
  
  #  with(as.list(c(state_initial,params)),{
  # manually extract disease states 
  
  M  = state_initial_age[[1]]
  
  S0 = state_initial_age[[2]]
  E0 = state_initial_age[[3]]
  A0 = state_initial_age[[4]]
  I0 = state_initial_age[[5]]
  R0 = state_initial_age[[6]]
  V0 = state_initial_age[[7]]
  
  S1 = state_initial_age[[8]]
  E1 = state_initial_age[[9]]
  A1 = state_initial_age[[10]]
  I1 = state_initial_age[[11]]
  R1 = state_initial_age[[12]]
  V1 = state_initial_age[[13]]
  
  S2 = state_initial_age[[14]]
  E2 = state_initial_age[[15]]
  A2 = state_initial_age[[16]]
  I2 = state_initial_age[[17]]
  R2 = state_initial_age[[18]]
  V2 = state_initial_age[[19]]
  
  S3 = state_initial_age[[20]]
  E3 = state_initial_age[[21]]
  A3 = state_initial_age[[22]]
  I3 = state_initial_age[[23]]
  R3 = state_initial_age[[24]]
  V3 = state_initial_age[[25]]
  
  Z  = state_initial_age[[26]]
  
  # extract parameters 
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
  age_r  = params[['age_r']]
  N      = sum(unlist(states))
  uptake = params[['uptake']]                     
  eff = params[['eff']] 
  waning_rate = params[['waning_rate']]

  # update age  of initial states 
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
  
  # get day of the year
  t1 = t %% 365
  
  # seasonal administration of intervention [end of October/start of RSV season, day=300]
  if (t1 == 300){
    dV0_vac = (S0 * uptake* eff)+(E0*uptake* eff)+(A0*uptake* eff)
    V0 = V0 + dV0_vac
    dV0 = dV0 + dV0_vac # update states

    dS0_vac = -S0*uptake * eff
    S0 = S0 + dS0_vac
    dS0 = dS0 + dS0_vac
    dE0_vac = -E0*uptake * eff
    E0 = E0 + dE0_vac
    dE0 = dE0 + dE0_vac
    dA0_vac = -A0*uptake * eff
    A0 = A0 + dA0_vac
    dA0 = dA0 + dA0_vac

    dV1_vac = (S1 * uptake* eff)+(E1*uptake* eff)+(A1*uptake* eff)
    V1 = V1 + dV1_vac
    dV1 = dV1 + dV1_vac

    dS1_vac = -S1*uptake * eff
    S1 = S1 + dS1_vac
    dS1 = dS1 + dS1_vac
    dE1_vac = -E1*uptake * eff
    E1 = E1 + dE1_vac
    dE1 = dE1 + dE1_vac
    dA1_vac = -A1*uptake * eff
    A1 = A1 + dA1_vac
    dA1 = dA1 + dA1_vac

    dV2_vac = (S2 * uptake* eff)+(E2*uptake* eff)+(A2*uptake* eff)
    V2 = V2 + dV2_vac
    dV2 = dV2 + dV2_vac

    dS2_vac = -S2*uptake * eff
    S2 = S2 + dS2_vac
    dS2 = dS2 + dS2_vac
    dE2_vac = -E2*uptake * eff
    E2 = E2 + dE2_vac
    dE2 = dE2 + dE2_vac
    dA2_vac = -A2*uptake * eff
    A2 = A2 + dA2_vac
    dA2 = dA2 + dA2_vac

    dV3_vac = (S3 * uptake* eff)+(E3 * uptake * eff)+(A3 * uptake * eff)
    V3 = V3 + dV3_vac
    dV3 = dV3 + dV3_vac

    dS3_vac = -S3*uptake * eff
    S3 = S3 + dS3_vac
    dS3 = dS3 + dS3_vac
    dE3_vac = -E3*uptake * eff
    E3 = E3 + dE3_vac
    dE3 = dE3 + dE3_vac
    dA3_vac = -A3*uptake * eff
    A3 = A3 + dA3_vac
    dA3 = dA3 + dA3_vac
  }

  # transmission rate for each day of year
  beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))
  
  # force of infection per exposure level (0,1,2,3)
  lambda0 <- beta_value * d0*(sum(I0)+(sum(A0)*alpha)) / N # d0=1, for uniformity, not in David's model 
  lambda1 <- beta_value * d1*(sum(I1)+(sum(A1)*alpha)) / N
  lambda2 <- beta_value * (d1*d2)*(sum(I2)+(sum(A2)*alpha)) / N
  lambda3 <- beta_value * (d1*d2*d3)*(sum(I3)+(sum(A3)*alpha)) / N
  
  # age adjusted ODEs
  dM  = dM + propR*mu - M*xi
  dS0 = dS0 + (1-propR)*mu + M*xi - S0*lambda0  + V0*waning_rate
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
  
  dZ = dZ + sigma*(E0+E1+E2+E3)
  
  return(list(c(dM,dS0,dE0,dA0,dI0,dR0,dV0,dS1,dE1,dA1,dI1,dR1,dV1,dS2,dE2,dA2,dI2,dR2,dV2,dS3,dE3,dA3,dI3,dR3,dV3,dZ)))
  #})
}
############################ sanity checks #####################
params <- list( 
  # t          = t,                                  # time horizon
  b1           = 1.997732,                            # mean seasonal transmission rate
  phi          = 0.2401791,                           # seasonal offset/peak transmission seasonal
  psi          = 0.615044,                            # width of seasonal peak/standard deviation
  sigma        = 1/4.863965,                          # susceptibility to primary infection
  d0           = 1,                                  # relative susceptibility to primary infection [d0=1 not included in Hodgson et al, added for uniformity]
  d1           = 0.8908731,                           # relative susceptibility to secondary infection, relative to  primary infection
  d2           = 0.7912774,                           # relative susceptibility to tertiary infection, relative to secondary infection
  d3           = 0.3428649,                           # relative susceptibility to subsequent infections after third infection, relative to tertiary infection
  alpha        = 0.6085629,                           # reduction in infectiousness of asymptomatic infections
  propR        = 0.3,                                # proportion of neonates born with protection
  mu           = 1/1863,                                # daily birth rate # 1/1863
  xi           = 1/131.2745,                          # rate of loss maternal immunity
  probA        = 0.5,                                # probability infection is asymptomatic
  g0           = 6.062045,                             # rate of loss of primary infectiousness
  g1           = 8.776286e-01 ,                           # proportional decrease between secondary and primary infection
  g2           = 0.8122186,                           # proportional decrease between tertiary and secondary infection
  gamma3       = 1,                                  # duration of infectiousness for all infections after tertiary infection
  omega        = 1/356.8321,                          # rate of loss of disease-induced/post-infection immunity
  age_r        = 1/(20*365),                              # aging rate
  waning_rate1 = 1/180,                              # intervention waning rate - age group 1, generic will be updated in scenario analysis
  waning_rate2 = 1/180,                              # intervention waning rate - age group 2, generic will be updated in scenario analysis
  waning_rate3 = 1/180,                              # intervention waning rate - age group 3, generic will be updated in scenario analysis
  uptake1      = 0.9,                                # uptake rate of intervention - age group 1, generic will be updated in scenario analysis
  uptake2      = 0.6,                                # uptake rate of intervention - age group 2, generic will be updated in scenario analysis
  uptake3      = 0.7,                                # uptake rate of intervention - age group 3, generic will be updated in scenario analysis
  eff1         = 0.78,                               # intervention efficacy - age group 1, generic will be updated in scenario analysis
  eff2         = 0.66,                               # intervention efficacy - age group 2, generic will be updated in scenario analysis
  eff3         = 0.65                                # intervention efficacy - age group 3, generic will be updated in scenario analysis
) 

state_initial_sanity <- c(M=c(100,50,20), S0=c(100,50,20), E0=c(100,50,20), A0=c(100,50,20), I0=c(100,50,20),  R0=c(100,50,20), V0 = c(0,0,0),
                       S1=c(100,50,20), E1=c(100,50,20), A1=c(100,50,20), I1=c(100,50,20),  R1=c(100,50,20), V1 = c(0,0,0),
                       S2=c(100,50,20), E2=c(100,50,20), A2=c(100,50,20), I2=c(100,50,20),  R2=c(100,50,20), V2 = c(0,0,0),
                       S3=c(100,50,20), E3=c(100,50,20), A3=c(100,50,20), I3=c(100,50,20),  R3=c(100,50,20), V3 = c(0,0,0),
                       Z=c(0,0,0))

# health = 0, aging_rate = 0, only change intervention (t=300)
params_check <- params
params_check[1:18] <- 0
params_check$age_r <- 0
params_check[20:28] <- c(0,0.4,1)
params_check$d0 <- 0
params_check$gamma0 <- 0
params_check$gamma1 <- 0
params_check$gamma2 <- 0
params_check$gamma3 <- 0
time_ode <- seq(0,365, 1)
matrix(unlist(rsv_model_age(300,state_initial_sanity,params_check)),ncol=3,byrow = T) 
run_check <- ode(y=state_initial_sanity, times = time_ode, func = rsv_model_age, parms = params_check)
head(run_check)
tail(run_check)

# age = 0, t≠300 (no intervention), only health changes 
params_check_1 <- params
params_check_1$age_r <- 0
state_initial_check = state_initial_sanity
state_initial_check[!grepl('..1', names(state_initial_check))] =0
matrix(unlist(rsv_model_age(1,state_initial_check,params_check_1)),ncol=3,byrow = T) 

run_check_1 <- ode(y=state_initial_age, times = time_ode, func = rsv_model_age, parms = params, method = 'ode23')
head(run_check_1)
tail(run_check_1)
# compare with initial states without age groups
state_initial <- c(M=170, S0=170, E0=170, A0=170, I0=170, R0=170, V0=0,
                          S1=170, E1=170, A1=170, I1=170, R1=170, V1=0,
                          S2=170, E2=170, A2=170, I2=170, R2=170, V2=0,
                          S3=170, E3=170, A3=170, I3=170, R3=170, V3=0,
                   Z=0)
run_check_original <- ode(y=state_initial, times = time_ode, func = rsv_model, parms = params, method = 'ode23')
rsv_model(1, state_initial, params)
head(run_check_original)
tail(run_check_original)

# t ≠ 300, aging_r = 0, health = 0, only waning 
params_check_2 <- params
params_check_2[1:18] <- 0
params_check_2[23:28] <- 0
params_check_2$d0 <- 0
params_check_2$age_r <- 0
params_check_2$gamma0 <- 0
params_check_2$gamma1 <- 0
params_check_2$gamma2 <- 0
params_check_2$gamma3 <- 0
params_check_2$waning_rate1 <- 1/180
params_check_2$waning_rate2 <- 1/80
params_check_2$waning_rate3 <- 1/20
state_initial_check <- state_initial_age
state_initial_check[grepl('V', names(state_initial_check))] = 100

matrix(unlist(rsv_model_age(1,state_initial_check,params_check_2)),ncol=3,byrow = T)
run_d <- ode(y=state_initial_check, times = time_ode, func = rsv_model_age, parms = params_check_2, method = 'ode23')
head(run_d[,"V01"])
# check issue with waning immune of vaxx

# run ode and generate plots 
library(deSolve)
# # 60 years 
# time_ode <- seq(0,(60*365), 1)
# run_ode <- ode(y=state_initial_age, times = time_ode, func = rsv_model_age, parms = params) 
# head(run_ode)
# # error 1: an excessive amount of work (> maxsteps ) was done, but integration was not successful - increase maxsteps
# # error 2: stiff ODEs
# unlist(params) # looks fine, possible issue with no. in initial states at t=0

# model run, no intevention 
time_ode <- seq(0,365, 1)
run_ode <- ode(y=current_state, times = time_ode, func = rsv_model_age, parms = params, method = 'ode23')
head(run_ode)
tail(run_ode)
plot(run_ode[, 8])
plot(run_ode[, "V21"]) 
plot(run_ode[, "V31"]) # looks to be some issue but there is movement 

sum()

boxplot(run_ode[, "Z1"], run_ode[, "Z2"], run_ode[, "Z3"])
hist(run_ode[, "Z2"])
hist(run_ode[, "Z3"])

plot(run_ode[, "Z1"])
plot(run_ode[, "Z2"])
plot(run_ode[, "Z3"])

hist(c(run_ode[, "Z1"],run_ode[, "Z2"], run_ode[, "Z3"]))


