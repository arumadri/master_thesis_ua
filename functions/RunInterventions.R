#### necessary functions ####

ParameterValuesforODE <- function(currentParamValues) {
  if(length(currentParamValues) != 25) {
    stop("currentParamValues must have 25 elements.")
  }
  
  parameterValuesTemp <- numeric(25)
  
  names(parameterValuesTemp) <- c("xi", "si", "ga0", "g1", "g2", "om", "pA1", "pA2", "pA3", "pA4", 
                                  "alpha_i", "d1", "d2", "d3", "phi", "qp", "qc", "b1", "psi", 
                                  "c5ep1", "c5ep2", "ep5", "ep6", "I1", "I2")
  
  parameterValuesTemp <- setNames(currentParamValues, names(parameterValuesTemp))
  
 # initialize 
  ep_t <- numeric(25)
  pA <- numeric(25)
  parameterValues <- parameterValuesTemp
  
  for (a in 1:16) {
    ep_t[a] <- exp(parameterValues["c5ep1"] + (a - 1) * parameterValues["c5ep2"])
  }
  for (a in 17:23) {
    ep_t[a] <- parameterValues["ep5"]
  }
  for (a in 24:25) {
    ep_t[a] <- parameterValues["ep6"]
  }
  
  pA[1:12] <- rep(parameterValues["pA1"], 12)
  pA[13:16] <- rep(parameterValues["pA2"], 4)
  pA[17:18] <- rep(parameterValues["pA3"], 2)
  pA[19:25] <- rep(parameterValues["pA4"], 7)
  
  assign("ep_t", ep_t, envir = .GlobalEnv)
  assign("pA", pA, envir = .GlobalEnv)
  assign("parameterValues", parameterValues, envir = .GlobalEnv)
  
}

initial_M <- function(parameterValues, ageStratification, populationPerAgeGroup) {
  xi <- 1.0 / parameterValues[['xi']]
  
  init_con <- numeric(0)
  
  for (i in 1:(length(ageStratification) - 1)) {
    
    cdf_lower <- pexp((365 * ageStratification[i]), rate = xi)
    cdf_upper <- pexp((365 * ageStratification[i + 1]), rate = xi)
    
    init_con_temp <- (cdf_upper - cdf_lower) / ((365 * ageStratification[i + 1] - 365 * ageStratification[i]) * xi)
    init_con[i] <- init_con_temp * populationPerAgeGroup[i]
  }
  
  cdf_last <- pexp(365 * 90, rate = xi) - pexp(365 * ageStratification[length(ageStratification)], rate = xi)
  last_age_group_init_con <- cdf_last / ((365 * 90 - 365 * ageStratification[length(ageStratification)]) * xi) * populationPerAgeGroup[length(ageStratification)]
  
  init_con <- c(init_con, last_age_group_init_con)
  
  return(init_con)
}

poisson_cdf <- function(l, a, x) {
  if (l == 0.0 || a == 0.0) {
    return(ppois(x, lambda = 0.000001))
  } else {
    return(ppois(x, lambda = l * a))
  }
}

initialProportionExposure <- function(l, a1, a2) {
  
  prop <- vector("numeric", 4)
  prop[1] <- abs(poisson_cdf(l, a2, 0) - poisson_cdf(l, a1, 0)) / ((a2 - a1) * l)
  prop[2] <- abs(poisson_cdf(l, a2, 1) - poisson_cdf(l, a1, 1)) / ((a2 - a1) * l)
  prop[3] <- abs(poisson_cdf(l, a2, 2) - poisson_cdf(l, a1, 2)) / ((a2 - a1) * l)
  prop[4] <- 1 - sum(prop[1:3])
  
  return(prop)
}

generateInitialStates <- function(parameterValues, ageStratification, populationPerAgeGroup) {
  
  populationMatPro <- initial_M(parameterValues, ageStratification, populationPerAgeGroup)
  
  I1 <- parameterValues[["l1"]]
  I2 <- parameterValues[["l2"]]
  I3 <- 0.5
  si <- 1.0 / parameterValues[["si"]]
  g0 <- 1.0 / parameterValues[["g0"]]
  g1 <- 1.0 / (parameterValues[["g0"]] * parameterValues[["g1"]])
  g2 <- 1.0 / (parameterValues[["g0"]] * parameterValues[["g1"]] * parameterValues[["g2"]])
  d1 <- parameterValues[["d1"]]
  d2 <- parameterValues[["d1"]] * parameterValues[["d2"]]
  d3 <- parameterValues[["d1"]] * parameterValues[["d2"]] * parameterValues[["d3"]]
  
  initialStates <- c()
  for (a in 1:25) {
    if (a < 25) {
      a1 <- ageStratification[a]
      a2 <- ageStratification[a + 1]
    } else {
      a1 <- ageStratification[a]
      a2 <- 90
    }
    propEachExposureGroup <- initialProportionExposure(I3, a1, a2)  
    pI1 <- propEachExposureGroup[1]
    pI2 <- propEachExposureGroup[2]
    pI3 <- propEachExposureGroup[3]
    pI4 <- propEachExposureGroup[4]
    age_size <- populationPerAgeGroup[a] - populationMatPro[a]
    
    initialStates_i <- c(populationMatPro[a],  # M group
                         
                         pI1 * age_size * (1 - I1) * (1 - I2),  # S0
                         pI1 * age_size * I1 * si / (si + g0),  # E0
                         pI1 * age_size * I1 * g0 / (si + g0) * pA[a],  # A0
                         pI1 * age_size * I1 * g0 / (si + g0) * (1 - pA[a]),  # I0
                         pI1 * age_size * (1 - I1) * I2,  # R0,
                         0, # V0
                         
                         pI2 * age_size * (1 - d1 * I1) * (1 - I2),  # S1
                         pI2 * age_size * d1 * I1 * si / (si + g1),  # E1
                         pI2 * age_size * d1 * I1 * g1 / (si + g1) * pA[a],  # A1
                         pI2 * age_size * d1 * I1 * g1 / (si + g1) * (1 - pA[a]),  # I1
                         pI2 * age_size * (1 - d1 * I1) * I2,  # R1
                         0, # V1
                         
                         pI3*age_size*(1 - d2*I1)*(1.0-I2),      # S2
                         pI3*age_size*d2*I1*si/(si+g2),            # E2
                         pI3*age_size*d2*I1*g2/(si+g2)*pA[a],      # A2
                         pI3*age_size*d2*I1*g2/(si+g2)*(1-pA[a]),  # I2
                         pI3*age_size*(1 - d2*I1)*I2,      # R2
                         0, # V2
                         
                         pI4*age_size*(1 - d3*I1)*(1-I2),        # S3
                         pI4*age_size*d3*I1*si/(si+g2),       # E3
                         pI4*age_size*d3*I1*g2/(si+g2)*pA[a], # A3
                         pI4*age_size*d3*I1*g2/(si+g2)*(1-pA[a]),   # I3
                         pI4*age_size*(1 - d3*I1)*I2, # R3
                         0,    # V3
                         
                         0  # Z
    )
    initialStates <- c(initialStates, initialStates_i) # these are vaccine states which are intially 0
  }
  return(initialStates)
}

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

plot_figures <- function(ode_output, programme_name, subtitle){
  
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
  }
  
  library(pacman)
  
  pacman::p_load(tidyverse, patchwork)
  
  ode_output <- as.data.frame(ode_output)
  programme_name <- as.character(programme_name)
  
  # p_maternal
  p_maternal <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("M1","M2", "M3"))  
  
  p_maternal <- ggplot(p_maternal, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Maternal",
         x = "Time",
         y = "Population") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,3650), breaks = seq(0,3650,365),
                       labels = c("0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"), 
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none") 
  # level 0
  p_level0 <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("S01", "S02", "S03", "E01", "E02", "E03", "A01", "A02", "A03", 
                        "I01", "I02", "I03", "R01", "R02", "R03", "V01", "V02", "V03")) %>%
    mutate(state = factor(state, levels = c("S01", "S02", "S03", "E01", "E02", "E03", "A01", "A02", "A03",
                                            "I01", "I02", "I03", "R01", "R02", "R03", "V01", "V02", "V03")))
  
  p_level0 <- ggplot(p_level0, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Exposure level 0",
         x ="Time (years)",
         y = "Population",
         tag = "A") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,3650), breaks = seq(0,3650,365),
                       labels = c("0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none") 
  
  # level 1
  p_level1 <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("S11", "S12", "S13", "E11", "E12", "E13", "A11", "A12", "A13", 
                        "I11", "I12", "I13", "R11", "R12", "R13", "V11", "V12", "V13"))%>%
    mutate(state = factor(state, levels = c("S11", "S12", "S13", "E11", "E12", "E13", "A11", "A12", "A13",
                                            "I11", "I12", "I13", "R11", "R12", "R13", "V11", "V12", "V13")))
  
  p_level1 <- ggplot(p_level1, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Exposure level 1",
         x = "Time (years)",
         y = "",
         tag = "B") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,3650), breaks = seq(0,3650,365),
                       labels = c("0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none") 
  
  # level 2
  p_level2 <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("S21", "S22", "S23", "E21", "E22", "E23", "A21", "A22", "A23", 
                        "I21", "I22", "I23", "R21", "R22", "R23", "V21", "V22", "V23"))%>%
    mutate(state = factor(state, levels = c("S21", "S22", "S23", "E21", "E22", "E23", "A21", "A22", "A23",
                                            "I21", "I22", "I23", "R21", "R22", "R23", "V21", "V22", "V23")))
  
  p_level2 <- ggplot(p_level2, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Exposure level 2",
         x = "Time (years)",
         y = "Population",
         tag = "C") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,3650), breaks = seq(0,3650,365),
                       labels = c("0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          plot.tag=element_text(face="bold"),
          legend.position = "none") 
  
  # level 3
  p_level3 <- ode_output %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("S31", "S32", "S33", "E31", "E32", "E33", "A31", "A32", "A33", 
                        "I31", "I32", "I33", "R31", "R32", "R33", "V31", "V32", "V33"))%>%
    mutate(state = factor(state, levels = c("S31", "S32", "S33", "E31", "E32", "E33", "A31", "A32", "A33",
                                            "I31", "I32", "I33", "R31", "R32", "R33", "V31", "V32", "V33")))
  
  p_level3 <- ggplot(p_level3, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Exposure level 3",
         x = "Time (years)",
         y = "",
         tag = "D") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(0,3650), breaks = seq(0,3650,365),
                       labels = c("0","1", "2", "3", "4", "5", "6", "7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=7, color = "black"),
          axis.text.y=element_text(family="sans",size=7, color = "black"),
          legend.text=element_text(family="sans",size=8, color = "black"),
          legend.title=element_text(family="sans",size=8, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  
  # joint 
  plot_base <- (p_level0 | p_level1) / (p_level2 | p_level3)
  
  plot_base <- plot_base + plot_annotation(
    title = programme_name,
    subtitle = subtitle,
    theme = theme(
      plot.title = element_text(family="sans",size=15, color = "black"),
      plot.subtitle = element_text(family="sans",size=13, color = "black")+
        theme_classic()
    )
  )
  
  return(plot_base)
}

plot_incidence <- function(ode_output, programme_name, subtitle){
  
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
  }
  
  library(pacman)
  
  pacman::p_load(tidyverse, patchwork)

  ode_output <- as.data.frame(ode_output)
  programme_name <- as.character(programme_name)
  # level 0
  results_base_stable <- as.data.frame(ode_output[2190:3650,])
  results_base_0 <- results_base_stable
  annual_incidence_0 <- results_base_0  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E01", "E02", "E03"))%>%
    mutate(state = factor(state, levels = c("E01", "E02", "E03")))
  
  annual_incidence_0 <- ggplot(annual_incidence_0, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 9900), data = subset(annual_incidence_0, state == "E01"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 489), data = subset(annual_incidence_0, state == "E02"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 2.63), data = subset(annual_incidence_0, state == "E03"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 7500), data = subset(annual_incidence_0, state == "E01"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 360), data = subset(annual_incidence_0, state == "E02"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 1.93), data = subset(annual_incidence_0, state == "E03"), color = "black", linetype = "dashed") +
    labs(title = "Exposure level 0",
         x = "Time (years)",
         y = "Population",
         tag = "A") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2190,3650), breaks = seq(2190,3650, 365),
                       labels = c("6","7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=13, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1, size = 13, color = "black"),
          axis.title.y = element_text(size = 13, colour = "black"),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 9, face = "bold"),
          legend.position = "none" ) 
  annual_incidence_0
  
  # level 1
  results_base_1 <- as.data.frame(results_base_stable)
  annual_incidence_1 <- results_base_1  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E11", "E12", "E13"))%>%
    mutate(state = factor(state, levels = c("E11", "E12", "E13")))
  
  annual_incidence_1 <- ggplot(annual_incidence_1, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 8465), data = subset(annual_incidence_1, state == "E11"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 2375), data = subset(annual_incidence_1, state == "E12"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 41.7), data = subset(annual_incidence_1, state == "E13"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 5885), data = subset(annual_incidence_1, state == "E11"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 1605), data = subset(annual_incidence_1, state == "E12"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 27.8), data = subset(annual_incidence_1, state == "E13"), color = "black", linetype = "dashed") +
    labs(title = "Exposure level 1",
         x = "Time (years)",
         y = "",
         tag = "B") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2190,3650), breaks = seq(2190,3650, 365),
                       labels = c("6","7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=13, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1, size = 13, color = "black"),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 9, face = "bold"),
          legend.position = "none" ) 
  
  annual_incidence_1
  # level 2
  results_base_2 <- as.data.frame(results_base_stable)
  annual_incidence_2 <- results_base_2  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E21", "E22", "E23"))%>%
    mutate(state = factor(state, levels = c("E21", "E22", "E23")))
  
  annual_incidence_2 <- ggplot(annual_incidence_2, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 8050), data = subset(annual_incidence_2, state == "E21"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 4670), data = subset(annual_incidence_2, state == "E22"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 131), data = subset(annual_incidence_2, state == "E23"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 4230), data = subset(annual_incidence_2, state == "E21"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 2380), data = subset(annual_incidence_2, state == "E22"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 65), data = subset(annual_incidence_2, state == "E23"), color = "black", linetype = "dashed") +
    labs(title = "Exposure level 2",
         x = "Time (years)",
         y = "Population",
         tag = "C") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2190,3650), breaks = seq(2190,3650, 365),
                       labels = c("6","7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=13, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1, size = 13, color = "black"),
          axis.title.y = element_text(size = 13, colour = "black"),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 9, face = "bold"),
          legend.position = "none" ) 
  
  annual_incidence_2
  
  # level 3
  results_base_3 <- as.data.frame(results_base_stable)
  annual_incidence_3 <- results_base_3  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E31", "E32", "E33"))%>%
    mutate(state = factor(state, levels = c("E31", "E32", "E33")))
  
  annual_incidence_3 <- ggplot(annual_incidence_3, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 1), data = subset(annual_incidence_3, state == "E31"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 1), data = subset(annual_incidence_3, state == "E32"), color = "black", linetype = "dotdash") +
    geom_hline(aes(yintercept = 1), data = subset(annual_incidence_3, state == "E33"), color = "black", linetype = "dotdash") +
    labs(title = "Exposure level 3",
         x = "Time (years)",
         y = "",
         tag = "D") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2190,3650), breaks = seq(2190,3650, 365),
                       labels = c("6","7", "8", "9", "10") ) +
    theme(axis.text.x=element_text(family="sans",size=12, color = "black"),
          axis.text.y=element_text(family="sans",size=12, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1, size = 13, color = "black"),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 9, face = "bold"),
          legend.position = "none" ) 
  
  annual_incidence_3
  
  annual_incidence <- (annual_incidence_0 | annual_incidence_1)/ (annual_incidence_2 | annual_incidence_3)
  
  
  annual_incidence <- annual_incidence + plot_annotation(
    title = programme_name,
    subtitle = subtitle,
    theme = theme(
      plot.title = element_text(family="sans",size=16, color = "black", face = "bold"),
      plot.subtitle = element_text(family="sans",size=14, color = "black")+
        theme_classic()
    )
  )
  
  return(annual_incidence)
}

plot_incidence_annual <- function(ode_output, programme_name, subtitle){
  
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
  }
  
  library(pacman)
  
  pacman::p_load(tidyverse, patchwork)
  
  ode_output <- as.data.frame(ode_output)
  programme_name <- as.character(programme_name)
  # level 0
  results_stable <- as.data.frame(ode_output[2565:2940,])
  results_stable_0 <- results_stable
  compare_incidence_0 <- results_stable_0  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E01", "E02", "E03"))%>%
    mutate(state = factor(state, levels = c("E01", "E02", "E03")))
  
  compare_incidence_0 <- ggplot(compare_incidence_0, aes(x = time, y = value, color = state)) +
    geom_line() +
    labs(title = "Exposure level 0",
         x = "Time (months)",
         y = "Population",
         tag = "A") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2565,2940), breaks = seq(2565,2940,90),
                       labels = c("Oct","Jan", "Apr", "Jul", "")  ) +
    theme(axis.text.x=element_text(family="sans",size=10, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  compare_incidence_0
  
  # level 1
  results_stable_1 <- as.data.frame(results_stable)
  compare_incidence_1 <- results_stable_1  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E11", "E12", "E13"))%>%
    mutate(state = factor(state, levels = c("E11", "E12", "E13")))
  
  compare_incidence_1 <- ggplot(compare_incidence_1, aes(x = time, y = value, color = state)) +
    geom_line() +
    labs(title = "Exposure level 1",
         x = "Time (months)",
         y = "",
         tag = "B") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2565,2940), breaks = seq(2565,2940,90),
                       labels = c("Oct","Jan", "Apr", "Jul", "")  )  +
    theme(axis.text.x=element_text(family="sans",size=10, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  
  compare_incidence_1
  # level 2
  results_stable_2 <- as.data.frame(results_stable)
  compare_incidence_2 <- results_stable_2  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E21", "E22", "E23"))%>%
    mutate(state = factor(state, levels = c("E21", "E22", "E23")))
  
  compare_incidence_2 <- ggplot(compare_incidence_2, aes(x = time, y = value, color = state)) +
    geom_line() +
    labs(title = "Exposure level 2",
         x = "Time (months)",
         y = "Population",
         tag = "C") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2565,2940), breaks = seq(2565,2940,90),
                       labels = c("Oct","Jan", "Apr", "Jul", "")  ) +
    theme(axis.text.x=element_text(family="sans",size=10, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  
  compare_incidence_2
  
  # level 3
  results_stable_3 <- as.data.frame(results_stable)
  compare_incidence_3 <- results_stable_3  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E31", "E32", "E33"))%>%
    mutate(state = factor(state, levels = c("E31", "E32", "E33")))
  
  compare_incidence_3 <- ggplot(compare_incidence_3, aes(x = time, y = value, color = state)) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1)) + 
    labs(title = "Exposure level 3",
         x = "Time (months)",
         y = "",
         tag = "D") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2565,2940), breaks = seq(2565,2940,90),
                       labels = c("Oct","Jan", "Apr", "Jul", "")  )  +
    theme(axis.text.x=element_text(family="sans",size=10, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  
  compare_incidence_3
  
  compare_incidence <- (compare_incidence_0 | compare_incidence_1)/ (compare_incidence_2 | compare_incidence_3)
  
  
  compare_incidence <- compare_incidence + plot_annotation(
    title = programme_name,
    subtitle = subtitle,
    theme = theme(
      plot.title = element_text(family="sans",size=16, color = "black", face = "bold"),
      plot.subtitle = element_text(family="sans",size=14, color = "black")+
        theme_classic()
    )
  )
  
  return(compare_incidence)
}

new_cases <- function(ode_output){
  start_scenario <- (6*365 + 274) # 1st October
  end_scenario <- 2830 # end of season 
  
  relevant_data <- ode_output[start_scenario:end_scenario, ]
  
  daily_changes <- diff(relevant_data[,c("Z1", "Z2", "Z3")])
  daily_changes_df <- as.data.frame(daily_changes)
  
  new_cases <- colSums(daily_changes)
  
  return(new_cases)
}
