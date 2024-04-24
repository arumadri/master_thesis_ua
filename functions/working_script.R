
# data frame 
age_labels <- c("0-1mo", "1-2mo", "2-3mo", "3-4mo", "4-5mo", "5-6mo", "6-7mo", "7-8mo", 
                "8-9mo", "9-10mo", "10-11mo", "11mo-1yr", "1-2yr", "2-3yr", "3-4yr", "4-5yr", "5-9yr", "10-14yr", 
                "15-24yr", "25-34yr", "35-44yr", "45-54yr", "55-64yr", "65-74yr", "75yr+")

data_for_plot <- data.frame(AgeRange = age_labels, Incidence = incidences)

library(ggplot2)

# a <- c(5, 11, 17, 23)  # Indices from which to sum values
# symptomatic <- numeric(25)  # Create an empty vector to store the aggregated data
# 
# for (i in 1:25) {
#   result_vector <- unlist(results[[i]])  # Unlist the results for the current age group
#   # Sum the values at indices specified by 'a' and assign to symptomatic[i]
#   symptomatic[i] <- sum(result_vector[a])
# }
# 
# symptomatic


# Set the factor levels for AgeRange to ensure they appear in the correct order
data_for_plot$AgeRange <- factor(data_for_plot$AgeRange, levels = age_labels)

results
# store annual incidences per age group
incidences_age <- list()

for (i in 1:length(results)) {
  incidences_age[[i]] <- results[[i]][366,27]  # 27th position (Z= incidence, 366 = total incidence at last day of the year)
}

unlist(incidences_age)

# data frame 
age_labels_age <- c("0-1mo", "1-2mo", "2-3mo", "3-4mo", "4-5mo", "5-6mo", "6-7mo", "7-8mo", 
                    "8-9mo", "9-10mo", "10-11mo", "11mo-1yr", "1-2yr", "2-3yr", "3-4yr", "4-5yr", "5-9yr", "10-14yr", 
                    "15-24yr", "25-34yr", "35-44yr", "45-54yr", "55-64yr", "65-74yr", "75yr+")

data_for_plot <- data.frame(AgeRange = age_labels_age, Incidence = unlist(incidences_age))
data_for_plot$AgeRange <- factor(data_for_plot$AgeRange, levels = age_labels_age)
# barplot
incidenc_pdf <- ggplot(data_for_plot, aes(x = AgeRange, y = Incidence)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Annual Incidence of RSV infections by Age Group",
       x = "Age Range",
       y = "Incidence") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14), 
        plot.title = element_text(size = 16, hjust = 0.5)) +  
  scale_y_continuous(labels = scales::comma) 

# ggsave("incidence_rsv.pdf", plot = incidenc_pdf, device = "pdf", 
#        width = 16, height = 9, units = "in", dpi = 600)

################# lambda trials keeping them here ##############
rsv_model_age <- function(t, state_initial_age, params, contact_matrix) {
  # ... [code to extract state variables] ...
  
  # Assuming the total population N is the sum of all compartments for each age group
  N <- colSums(matrix(unlist(state_initial_age), nrow = length(state_initial_age) / 5, byrow = TRUE))
  
  # beta
  beta_value <- 1 + b1 * exp(-((t - phi)^2 / (2 * psi^2)))
  
  # Calculate the force of infection for each age group and infection level
  lambda <- matrix(0, nrow = 4, ncol = 5)  # 4 infection levels by number of age groups
  
  # Infectious individuals by age group and infection level
  I <- matrix(c(I0, I1, I2, I3), nrow = 4, byrow = TRUE)  # Each row corresponds to a different infection level
  A <- matrix(c(A0, A1, A2, A3), nrow = 4, byrow = TRUE)  # Each row corresponds to a different infection level
  
  # Calculate the force of infection
  for (level in 1:nrow(lambda)) {
    for (age_group in 1:ncol(lambda)) {
      infectious_pressure <- sum(contact_matrix[age_group,] * (I[level,] + alpha * A[level,]))
      lambda[level, age_group] <- beta_value * infectious_pressure / N[age_group]
    }
  }
  
  # ... [rest of the model code, using lambda for each age group and infection level] ...
  
  # Return the rate of change of state variables
}


rsv_model_age <- function(t, state_initial_age, params, contactMatrixPhy, contactMatrixCon) {
  
  beta_value <- (1 + b1 * exp(-((t - phi)^2 / (2 * psi^2))))
  
  # Number of asymptomatic and symptomatic individuals in each age group
  A <- c(sum(A0), sum(A1), sum(A2), sum(A3), sum(A4))
  I <- c(sum(I0), sum(I1), sum(I2), sum(I3), sum(I4))
  
  # Calculate the force of infection for each age group
  lambda <- matrix(0, nrow = 4, ncol = 5)  # Initialize lambda for each age group
  
  for (a in 1:length(N)) {
    contact_sum <- 0
    # Summation over all age groups
    for (b in 1:length(N)) {
      # Contact rate for physical and conversational between age groups a and b
      p_ab <- contactMatrixPhy[[b]][a]  # Adjust indexing depending on the structure of contact matrices
      c_ab <- contactMatrixCon[[b]][a]  # Adjust indexing depending on the structure of contact matrices
      
      contact_sum <- contact_sum + ((p_ab * qp + c_ab * qc) * (A[b] * alpha + I[b]) / N[b])
    }
    
    # Incorporate seasonal forcing into the force of infection for age group a
    lambda[a] <- seasonal_forcing * contact_sum
  }
  
  # ... [rest of the model code, using lambda for each age group] ...
  
  # Return the rate of change of state variables
}


beta_value <- (1 + b1 * exp(-((t - phi)^2 / (2 * psi^2))))
lambda <- matrix(0, nrow = 4, ncol = 5)  # Initialize lambda for each age group
delta_i <- c(d0,d1,d2,d3)

for (d in 1:delta_i) { # infection levels
  # Number of asymptomatic and symptomatic individuals in each age group for this level
  # Infectious individuals by age group and infection level
  I <- matrix(c(I0, I1, I2, I3), nrow = 4, byrow = TRUE)  # Each row corresponds to a different infection level
  A <- matrix(c(A0, A1, A2, A3), nrow = 4, byrow = TRUE)  # Each row corresponds to a different infection level
  
  for (a in 1:length(N)) {
    contact_sum <- 0
    # Summation over all age groups
    for (b in 1:length(N)) {
      # Contact rate for physical and conversational between age groups a and b
      p_ab <- contactMatrixPhy[[b]][a]  # Adjust indexing depending on the structure of contact matrices
      c_ab <- contactMatrixCon[[b]][a]  # Adjust indexing depending on the structure of contact matrices
      
      contact_sum <- contact_sum + ((p_ab * qp + c_ab * qc) * (A[b] * alpha + I[b]) / N[b])
    }
    
    # Incorporate seasonal forcing and relative susceptibility into the force of infection
    lambda[level, a] <- seasonal_forcing * contact_sum * delta_i[d]
  }
}
#################
# beta
beta_value = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))))

# lambda
lambda0 <- beta_value * d0 *(sum(I0)+(sum(A0)*alpha)) / N
lambda1 <- beta_value * d1* (sum(I1)+(sum(A1)*alpha)) / N
lambda2 <- beta_value * (d1*d2)*(sum(I2)+(sum(A2)*alpha)) / N
lambda3 <- beta_value * (d1*d2*d3)*(sum(I3)+(sum(A3)*alpha)) / N

######

# scenarios
# try 1 year
############################ run scenarios #####################
# source (Hodgson et al 2020 and 2024)

# base (pmab) scenario [Palivizimab for infants]
params_pmab <- params
params_pmab$waning_rate1 <- 1/30                                  
params_pmab$waning_rate2 <- 0                                  
params_pmab$waning_rate3 <- 0                                  
params_pmab$uptake1      <- 0.9                               
params_pmab$uptake2      <- 0                                
params_pmab$uptake3      <- 0                                
params_pmab$eff1         <- 0.338                                  
params_pmab$eff2         <- 0                                    
params_pmab$eff3         <- 0

# run ode and plot 
run_ode_pmab <- ode(y=state_initial_age, times = time_ode, func = rsv_model_age, parms = params, method = 'ode23') 
head(run_ode_pmab)
tail(run_ode_pmab)
plot(run_ode_pmab[, "Z1"])
plot(run_ode_pmab[, "Z2"]) 
plot(run_ode_pmab[, "Z3"]) # looks to be some issue but there is movement 


# nmab scenario [Niservimab for infants]
params_nmab <- params
params_nmab$waning_rate1 <- 1/150                                  
params_nmab$waning_rate2 <- 0                                  
params_nmab$waning_rate3 <- 0                                  
params_nmab$uptake1      <- 0.9*                             
params_nmab$uptake2      <- 0                                
params_nmab$uptake3      <- 0                                
params_nmab$eff1         <- 0.781                                   
params_nmab$eff2         <- 0                                    
params_nmab$eff3         <- 0

# mat scenario [Maternal immunization]
params_mat <- params
params_mat$waning_rate1 <- 1/180                                 
params_mat$waning_rate2 <- 1/180                                 
params_mat$waning_rate3 <- 0                                  
params_mat$uptake1      <- 0.6                               
params_mat$uptake2      <- 0.6*                                
params_mat$uptake3      <- 0                                
params_mat$eff1         <- 0.492                                  
params_mat$eff2         <- 0.659                                   
params_mat$eff3         <- 0

# old scenario [Older adult vaccination]
params_old <- params
params_old$waning_rate1 <- 0                                  
params_old$waning_rate2 <- 0                                
params_old$waning_rate3 <- 1/180                                  
params_old$uptake1      <- 0                                
params_old$uptake2      <- 0                               
params_old$uptake3      <- 0.7                                
params_old$eff1         <- 0                                   
params_old$eff2         <- 0                                    
params_old$eff3         <- 0.667

# mat_old scenario [maternal & older adult vaccination]
params_mat_old <- params
params_mat_old$waning_rate1 <- 1/180                                  
params_mat_old$waning_rate2 <- 1/180                                  
params_mat_old$waning_rate3 <- 1/180                                 
params_mat_old$uptake1      <- 0.6*                             
params_mat_old$uptake2      <- 0.6*                             
params_mat_old$uptake3      <- 0.7                                
params_mat_old$eff1         <- 0.492                                
params_mat_old$eff2         <- 0.659                                   
params_mat_old$eff3         <- 0.667
