############################ initial states ###################
rm(list = ls())
state_initial <- c(M=10e6, S0=20e6,  E0=5e6,   A0=1.5e6,  I0=3.5e6, R0=0,
                           S1=10e6,  E1=2.5e6, A1=7.5e5,  I1=1.3e6, R1=0,
                           S2=5e6,   E2=1.3e6, A2=3.8e5,  I2=6e5,   R2=0,
                           S3=2.5e6, E3=6e5,   A3=1.28e5, I3=3e5,   R3=0,
                   Z=0)
############################ parameters #####################
params <- list(b1         = 1.99626,                            # mean seasonal transmission rate
                   phi    = 0.621604,                           # seasonal offset/peak transmission seasonal
                   psi    = 0.25852,                            # width of seasonal peak/standard deviation
                   sigma  = 1/4.20082,                          # susceptibility to primary infection 
                   d1     = 0.881906,                           # relative susceptibility to secondary infection, relative to  primary infection
                   d2     = 0.622153,                           # relative susceptibility to tertiary infection, relative to secondary infection
                   d3     = 0.437988,                           # relative susceptibility to subsequent infections after third infection, relative to tertiary infection
                   alpha  = 0.486876,                           # reduction in infectiousness of asymptomatic infections
                   propR  = 0.3,                                # proportion of neonates born with protection
                   mu     = 1/1863,                             # daily birth rate
                   xi     = 1/123.886,                          # rate of loss maternal immunity
                   probA  = 0.5,                                # probability infection is asymptomatic
                   g0     = 5.6174,                             # rate of loss of primary infectiousness
                   g1     = 0.895885,                           # proportional decrease between secondary and primary infection
                   g2     = 0.910452,                           # proportional decrease between tertiary and secondary infection
                   gamma3 = 1,                                  # duration of infectiousness for all infections after tertiary infection
                   omega  = 1/353.187,                          # rate of loss of post-infection immunity
                   # l1     = 0.315599,                         # proportion of persons infected 
                   # l2     = 0.376428,                         # proportion of persons uninfected and protected
                   # qp     = 0.0970337,                        # probability of transmission per physical contact 
                   # qc     = 0.999361,                         # probability of transmission per conversational contact
                   age_r  = 1/(20*365)                          # aging rate
)
############################ function to calculate gamma #####################
calculate_gamma <- function(params) {
  
  # parameters from list
  g0 = params[['g0']]
  g1 = params[['g1']]
  g2 = params[['g2']]
  
  # new parameters
  params[['gamma0']] <<- 1/g0
  params[['gamma1']] <<- 1/(g0*g1)
  params[['gamma2']] <<- 1/(g0*g1*g2)
  
}

calculate_gamma(params)
############################ rsv model function #####################
rsv_model <- function(t,state_initial, params){
  # explicitly define states
  M  = state_initial['M']
  
  S0 = state_initial['S0']
  E0 = state_initial['E0']
  A0 = state_initial['A0']
  I0 = state_initial['I0']
  R0 = state_initial['R0']
  
  S1 = state_initial['S1']
  E1 = state_initial['E1']
  A1 = state_initial['A1']
  I1 = state_initial['I1']
  R1 = state_initial['R1']
  
  S2 = state_initial['S2']
  E2 = state_initial['E2']
  A2 = state_initial['A2']
  I2 = state_initial['I2']
  R2 = state_initial['R2']
  
  S3 = state_initial['S3']
  E3 = state_initial['E3']
  A3 = state_initial['A3']
  I3 = state_initial['I3']
  R3 = state_initial['R3']
  
  Z  = state_initial['Z']
  
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
  d1     = params[['d1']]
  d2     = params[['d2']]
  d3     = params[['d3']]
  N      = sum(unlist(state_initial))
  
  t1 = t %% 365
  beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))
  
  # lambda
  lambda0 <- beta_value * (I0+(A0*alpha)) / N
  lambda1 <- beta_value * d1* (I1+(A1*alpha)) / N
  lambda2 <- beta_value * (d1*d2)*(I2+(A2*alpha)) / N
  lambda3 <- beta_value * (d1*d2*d3)*(I3+(A3*alpha)) / N
  
  # ODEs
  dM  = - M*xi + propR*mu 
  
  dS0 = (1-propR)*mu + M*xi- S0*lambda0  
  dE0 = S0*lambda0 - E0*sigma 
  dA0 = E0*probA*sigma - A0*gamma0   
  dI0 = E0*(1-probA)*sigma - I0*gamma0 
  dR0 = I0*gamma0 + A0*gamma0 - R0*omega 
  
  dS1 = R0*omega - S1*lambda1 
  dE1 = S1*lambda1 - E1*sigma
  dA1 = E1*probA*sigma - A1*gamma1
  dI1 = E1*(1-probA)*sigma - I1*gamma1
  dR1 = I1*gamma1 + A1*gamma1 - R1*omega
  
  dS2 = R1*omega - S2*lambda2
  dE2 = S2*lambda2 - E2*sigma
  dA2 = E2*probA*sigma - A2*gamma2
  dI2 = E2*(1-probA)*sigma - I2*gamma2
  dR2 = I2*gamma2 + A2*gamma2 - R2*omega
  
  dS3 = R2*omega + R3*omega - S3*lambda3 
  dE3 = S3*lambda3 - E3*sigma
  dA3 = E3*probA*sigma - A3*gamma3
  dI3 = E3*(1-probA)*sigma - I3*gamma3
  dR3 = I3*gamma3 + A3*gamma3 - R3*omega
  
  dZ = sigma*(E0+E1+E2+E3)
  
  return(list(c(dM,dS0,dE0,dA0,dI0,dR0,dS1,dE1,dA1,dI1,dR1,dS2,dE2,dA2,dI2,dR2,dS3,dE3,dA3,dI3,dR3,dZ)))
}

rsv_model(1,state_initial,params)
state_initial
state_initial + unlist(rsv_model(1,state_initial,params))
# AGE 
############################ add age adjustment ##################### 
state_initial_age <- c(M=c(1,0,3), S0=c(2,1,1), E0=c(3,1,5), A0=c(1,2,1), I0=c(3,1,2),  R0=c(1,2,1),
                                   S1=c(2,1,1), E1=c(3,1,5), A1=c(1,2,1), I1=c(3,1,2),  R1=c(1,2,1),
                                   S2=c(2,1,1), E2=c(3,1,5), A2=c(1,2,1), I2=c(3,1,2),  R2=c(1,2,1),
                                   S3=c(2,1,1), E3=c(3,1,5), A3=c(1,2,1), I3=c(3,1,2),  R3=c(1,2,1),
                       Z=c(1,0,3))
state_initial_age
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
  
  S1 = state_initial_age[grepl('S1',names(state_initial_age))]
  E1 = state_initial_age[grepl('E1',names(state_initial_age))]
  A1 = state_initial_age[grepl('A1',names(state_initial_age))]
  I1 = state_initial_age[grepl('I1',names(state_initial_age))]
  R1 = state_initial_age[grepl('R1',names(state_initial_age))]
  
  S2 = state_initial_age[grepl('S2',names(state_initial_age))]
  E2 = state_initial_age[grepl('E2',names(state_initial_age))]
  A2 = state_initial_age[grepl('A2',names(state_initial_age))]
  I2 = state_initial_age[grepl('I2',names(state_initial_age))]
  R2 = state_initial_age[grepl('R2',names(state_initial_age))]
  
  S3 = state_initial_age[grepl('S3',names(state_initial_age))]
  E3 = state_initial_age[grepl('E3',names(state_initial_age))]
  A3 = state_initial_age[grepl('A3',names(state_initial_age))]
  I3 = state_initial_age[grepl('I3',names(state_initial_age))]
  R3 = state_initial_age[grepl('R3',names(state_initial_age))]
  
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
  d1     = params[['d1']]
  d2     = params[['d2']]
  d3     = params[['d3']]
  age_r  = params[['age_r']]
  N      = sum(unlist(state_initial_age))
  
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
  
  ageZ = update_ages(Z,age_r)
  dZ   = ageZ$dP
  Z    = ageZ$P
  
  # beta
  t1 = t %% 365
  beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))
  
  # lambda
  lambda0 <- beta_value * (sum(I0)+(sum(A0)*alpha)) / N
  lambda1 <- beta_value * d1* (sum(I1)+(sum(A1)*alpha)) / N
  lambda2 <- beta_value * (d1*d2)*(sum(I2)+(sum(A2)*alpha)) / N
  lambda3 <- beta_value * (d1*d2*d3)*(sum(I3)+(sum(A3)*alpha)) / N
  
  # age adjusted health - ODEs
  dM  = propR*mu*c(1,0,0) - M*xi
  
  dS0 = (1-propR)*mu*c(1,0,0) + M*xi - S0*lambda0  
  dE0 = S0*lambda0 - E0*sigma  
  dA0 = E0*probA*sigma - A0*gamma0 
  dI0 = E0*(1-probA)*sigma - I0*gamma0 
  dR0 = I0*gamma0 + A0*gamma0 - R0*omega 
  
  dS1 = R0*omega - S1*lambda1 
  dE1 = S1*lambda1 - E1*sigma
  dA1 = E1*probA*sigma - A1*gamma1
  dI1 = E1*(1-probA)*sigma - I1*gamma1
  dR1 = I1*gamma1 + A1*gamma1 - R1*omega
  
  dS2 = R1*omega - S2*lambda2
  dE2 = S2*lambda2 - E2*sigma
  dA2 = E2*probA*sigma - A2*gamma2
  dI2 = E2*(1-probA)*sigma - I2*gamma2
  dR2 = I2*gamma2 + A2*gamma2 - R2*omega
  
  dS3 = R2*omega + R3*omega - S3*lambda3 
  dE3 = S3*lambda3 - E3*sigma
  dA3 = E3*probA*sigma - A3*gamma3
  dI3 = E3*(1-probA)*sigma - I3*gamma3
  dR3 = I3*gamma3 + A3*gamma3 - R3*omega
  
  dZ = sigma*(E0+E1+E2+E3)
  
  return(list(c(dM,dS0,dE0,dA0,dI0,dR0,dS1,dE1,dA1,dI1,dR1,dS2,dE2,dA2,dI2,dR2,dS3,dE3,dA3,dI3,dR3,dZ)))
  #})
}
# model run 
rsv_model_age(1,state_initial_age,params)

state_initial_age
state_initial_age + unlist(rsv_model_age(1,state_initial_age,params))

# sanity check 
params_check <- params
params_check$age_r <- 0
state_initial_check <- state_initial_age*0
pattern <- "M1|S01|E01|I01|A01|R01|S11|E11|I11|A11|R11|S21|E21|I21|A21|R21|S31|E31|I31|A31|R31|Z1"
state_initial_check[grepl(pattern, names(state_initial_check))] <- state_initial
state_initial_check

# compare 
state_initial + unlist(rsv_model(1,state_initial,params))
state_initial_check + unlist(rsv_model_age(1,state_initial_check,params_check))
