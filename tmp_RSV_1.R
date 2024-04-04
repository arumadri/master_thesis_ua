rm(list = ls())
# define states
state_initial <- c(M =1, S=3,I = 4)

# define function
sir_model <- function(state_initial){
  # with(as.list(c(state_initial)),{
  # explicitly define M, S, I
  M = state_initial['M']
  S = state_initial['S']
  I = state_initial['I']

  
  # adjust health states
  dM = - M*0.01
  dS = M*0.01 - 0.03*I*S
  dI = 0.03*I*S
  
  # return
  return(list(c(dM,dS,dI)))
  # })
}

# call function
sir_model(state_initial)
state_initial
state_initial + unlist(sir_model(state_initial))


## AGE

state_initial_age <- c(M=c(1,0,3), S=c(2,1,1),I = c(3,1,5))
state_initial_age

update_ages <- function(pop,age_r){
  age_pop = pop*age_r
  dP = - age_pop +  c(0,age_pop[-length(age_pop)])
  P = pop + dP
  return(list(P=P,dP=dP))
}

# define function
sir_model_age <- function(state_intial_age,params){
#  with(as.list(c(state_initial,params)),{
  M = state_initial_age[grepl('M',names(state_initial_age))]
  S = state_initial_age[grepl('S',names(state_initial_age))]
  I = state_initial_age[grepl('I',names(state_initial_age))]
  
  age_r <- params['age_r']
  
  # age
  ageM = update_ages(M,age_r)
  dM   = ageM$dP
  M    = ageM$P
  
  ageS = update_ages(S,age_r)
  dS   = ageS$dP
  S    = ageS$P
  
  ageI = update_ages(I,age_r)
  dI   = ageI$dP
  I    = ageI$P
  
  # health
  dM = dM - M*0.01
  dS = dS + M*0.01 - 0.03*I*S
  dI = dI + 0.03*I*S
  
  return(list(c(dM,dS,dI)))
  #})
}

params <- c(age_r = 1/(20*365))
sir_model_age(state_initial_age,params)
  
state_initial_age
state_initial_age + unlist(sir_model_age(state_initial_age,params))

state_initial_age
state_initial_age + unlist(sir_model_age(state_initial_age,params))

# sanity check
params_check <- params * 0
state_initial_age <- state_initial_age*0
state_initial_age[grepl(1,names(state_initial_age))] <- state_initial
state_initial_age

# compare 
state_initial + unlist(sir_model(state_initial))
state_initial_age + unlist(sir_model_age(state_initial_age,params_check))

############################ beta overtime #####################
# beta_values <- numeric(2*365)
# for (t in 1:(2 * 365)) {
#   b1     = parameters[['b1']]
#   phi    = parameters[['phi']]
#   psi    = parameters[['psi']]
#   
#   t1 = t %% 365
#   beta_values[t] = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))
# }
# plot(beta_values, type = 'l')

# model overtime 
# for (t in 1:(2 * 365)) {
#   daily_beta <- beta_values[t]
#   
#   new_state <- rsv_model(t, state_initial, daily_beta)
# }


