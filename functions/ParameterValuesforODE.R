####
ParameterValuesforODE <- function(currentParamValues) {
  if(length(currentParamValues) != 25) {
    stop("currentParamValues must have 25 elements.")
  }
  
  parameterValuesTemp <- numeric(25)
  
  names(parameterValuesTemp) <- c("xi", "si", "ga0", "g1", "g2", "om", "pA1", "pA2", "pA3", "pA4", 
                                  "alpha_i", "d1", "d2", "d3", "phi", "qp", "qc", "b1", "psi", 
                                  "c5ep1", "c5ep2", "ep5", "ep6", "I1", "I2")
  
  parameterValuesTemp <- setNames(currentParamValues, names(parameterValuesTemp))
  
  # Define and compute ep_t and pA
  ep_t <- numeric(25)
  pA <- numeric(25)
  
  for (a in 1:16) {
    ep_t[a] <- exp(parameterValuesTemp["c5ep1"] + (a - 1) * parameterValuesTemp["c5ep2"])
  }
  for (a in 17:23) {
    ep_t[a] <- parameterValuesTemp["ep5"]
  }
  for (a in 24:25) {
    ep_t[a] <- parameterValuesTemp["ep6"]
  }
  
  pA[1:12] <- rep(parameterValuesTemp["pA1"], 12)
  pA[13:16] <- rep(parameterValuesTemp["pA2"], 4)
  pA[17:18] <- rep(parameterValuesTemp["pA3"], 2)
  pA[19:25] <- rep(parameterValuesTemp["pA4"], 7)
  
  list(ep_t = ep_t, pA = pA, parameterValues = parameterValuesTemp)
}
