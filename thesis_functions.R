# necessary package
install.packages("deSolve")
library(deSolve)
# ########################################## solve ODE ###############################################
dynamic_rsv_model <- function(propR = propR, # proportion of neonates born with protection
                              mu = mu,       # daily birth rate
                              xi = xi,       # rate of loss maternal immunity
                              beta = beta,   # force of infection
                              prob = prob,   # probability RSV infection is asymptomatic
                              delta = delta, # rate of loss exposure to infection
                              gamma = gamma, # rate of loss of infectiousness
                              omega = omega, # rate of loss of post-infection immunity
                              l1 = l1,       #
                              l2 = l2        #
                              ){
  # run M-SEAIR ODEs
  dynamic_m_seair <- function(times, states, parameters){
    with(as.list(c(states, parameters)), {
      N = M+S+E+I+R
      dM = propR*mu*N - xi*M
      dS = (1-propR)*mu*N + xi*M - beta*S + omega*R
      dE = beta*S - delta*E
      dA = prob*delta*E - gamma*A
      dI = (1-prob)*delta*E - gamma*I
      dR = gamma*A + gamma*I - omega*R
      return(list(c(dM, dS, dE, dA, dI, dR))) 
    })
  }
}
# initial values
start <- c(N=100,000,
           M=N*propR*xi, 
           S=N*(1-propR*xi)*propR*(1-l1)*(1-l2), 
           E=N*(1-propR*xi)*propR*(delta/(gamma+delta))*l1, 
           A=N*(1-propR*xi)*propR*(gamma/(gamma+delta))*propR*l1,
           I=N*(1-propR*xi)*propR*(gamma/(gamma+delta))*(1-propR)*l1,
           R=N*(1-propR*xi)*propR*(1-l1)*l2)
# parameters
parameters <- c(propR=1, mu=274, xi=1/60, beta=, prob=, delta=1/4.98, gamma=1/6.16, omega=1/358.9, l1=, l2=)
# time steps
times <- seq(1,365,1)
# model run 
run_ode <- ode(y = start, times = times, parms = parameters, func = dynamic_m_seair)
