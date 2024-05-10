sir_vaccination_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Implement vaccination starting at day 300
    t1 = time %% 365
    
    if (time == 300) {
      vaccinated_per_day <- nu * S
    } else {
      vaccinated_per_day <- 0
    }
    
    # The rate at which susceptibles get infected per day.
    infection_rate <- beta * I * S / N
    
    # Differential equations
    dS <- -infection_rate - vaccinated_per_day
    dV <- vaccinated_per_day - beta_v * V * I / N  # Assuming vaccinated individuals can still get infected but at a lower rate
    dI <- infection_rate - gamma * I
    dR <- gamma * I
    
    # Return the rate of change
    list(c(dS, dV, dI, dR))
  })
}

initial_state <- c(S = 990, V = 0, I = 10, R = 0)  # 990 susceptible, 0 vaccinated, 10 infected, 0 recovered
parameters <- c(beta = 0.3, gamma = 0.1, nu = 0.05, beta_v = 0.05)  # Transmission rate, recovery rate, vaccination rate, reduced transmission rate for vaccinated

# Total population
parameters["N"] <- sum(initial_state)

times <- seq(0, 5 * 365, by = 1)  # Time from day 0 to day 5*365

# Solve the ODEs
run_ode <- ode(y = initial_state, times = times, func = sir_vaccination_model, parms = parameters, method = 'ode23')

plot(run_ode[, 'V'])
