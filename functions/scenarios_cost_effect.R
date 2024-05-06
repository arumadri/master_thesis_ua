# scenarios
# try 1 year
############################ run scenarios #####################
# source (Hodgson et al 2020 and 2024)
library(deSolve)

time_ode <- seq(0, 3650, 1) 

# base: no intervention
params_none <- params
params_none$waning_rate1 <- 0                                  
params_none$waning_rate2 <- 0                                  
params_none$waning_rate3 <- 0                                  
params_none$uptake1      <- 0                              
params_none$uptake2      <- 0                                
params_none$uptake3      <- 0                                
params_none$vac_eff1         <- 0                                 
params_none$vac_eff2         <- 0                                    
params_none$vac_eff3         <- 0

results_base <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_none, method = "ode23")
# no. of cases # c()
round(new_cases(results_base), digits = 0) # c(7687394, 2505202, 53501)

# averted c(0,0,0)

#plot
model_stabilize <- plot_figures(results_base, programme_name = "Base: No intervention")
ggsave("fig_1.pdf", plot = model_stabilize, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

# pmab scenario [Palivizimab for infants]
params_pmab <- params
params_pmab$waning_rate1 <- 1/30                                  
params_pmab$waning_rate2 <- 0                                  
params_pmab$waning_rate3 <- 0                                  
params_pmab$uptake1      <- 0.9                              
params_pmab$uptake2      <- 0                                
params_pmab$uptake3      <- 0                                
params_pmab$vac_eff1         <- 0.338                                  
params_pmab$vac_eff2         <- 0                                    
params_pmab$vac_eff3         <- 0

results_pmab <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_pmab, method = "ode23")
# no. of cases 
round(new_cases(results_pmab), digits = 0) # c(7703600, 2508361, 53578)

# averted c(16206, 3159, 77)


#plot
plot_figures(results_pmab, programme_name = "Palivizimab for infants")

# nmab scenario [Niservimab for infants]
params_nmab <- params
params_nmab$waning_rate1 <- 1/150                                  
params_nmab$waning_rate2 <- 0                                  
params_nmab$waning_rate3 <- 0                                  
params_nmab$uptake1      <- 0.9                            
params_nmab$uptake2      <- 0                                
params_nmab$uptake3      <- 0                                
params_nmab$vac_eff1         <- 0.781                                   
params_nmab$vac_eff2         <- 0                                    
params_nmab$vac_eff3         <- 0

results_nmab <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_nmab, method = "ode23")
# no. of cases 
round(new_cases(results_nmab), digits = 0) # c(7720458, 2512249, 53820)

# averted c(33064, 7047, 319)

#plot
plot_figures(results_nmab, programme_name = "Niservimab for infants")

# mat scenario [Maternal immunization]
params_mat <- params
params_mat$waning_rate1 <- 1/180                                 
params_mat$waning_rate2 <- 1/180                                 
params_mat$waning_rate3 <- 0                                  
params_mat$uptake1      <- 0.6*0.5                            
params_mat$uptake2      <- 0.6*0.5                             
params_mat$uptake3      <- 0                                
params_mat$vac_eff1         <- 0.513                                 
params_mat$vac_eff2         <- 0.717                                   
params_mat$vac_eff3         <- 0

results_mat <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_mat, method = "ode23")

# no. of cases 
round(new_cases(results_mat), digits = 0) # c(7699588, 2552178, 53895)

# averted c(12194,46976,394)

#plot
plot_figures(results_mat, programme_name = "Maternal vaccine for mothers and infants")

# old scenario [Older adult vaccination]
params_old <- params
params_old$waning_rate1 <- 0                                  
params_old$waning_rate2 <- 0                                
params_old$waning_rate3 <- 1/180                                  
params_old$uptake1      <- 0                                
params_old$uptake2      <- 0                               
params_old$uptake3      <- 0.7                                
params_old$vac_eff1         <- 0                                   
params_old$vac_eff2         <- 0                                    
params_old$vac_eff3         <- 0.667

results_old <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                   parms = params_old, method = "ode23")
# no. of cases 
round(new_cases(results_old), digits = 0) # c(7687571, 2505379, 53978)

# averted c(177, 177, 477)

#plot
plot_figures(results_old, programme_name = "Elderly adult vaccination")

# # mat_old scenario [maternal & older adult vaccination]
# params_mat_old <- params
# params_mat_old$waning_rate1 <- 1/180                                  
# params_mat_old$waning_rate2 <- 1/180                                  
# params_mat_old$waning_rate3 <- 1/180                                 
# params_mat_old$uptake1      <- 0.6*0.5                             
# params_mat_old$uptake2      <- 0.6*0.5                           
# params_mat_old$uptake3      <- 0.7                                
# params_mat_old$vac_eff1         <- 0.513                                
# params_mat_old$vac_eff2         <- 0.717                                   
# params_mat_old$vac_eff3         <- 0.667
# 
# results_mat_old <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
#                    parms = params_mat_old, method = "ode23")
# # no. of cases averted
# round(new_cases(results_mat_old), digits = 0) # c(7699849, 2552431, 55344)
# 
# # averted c(12455, 47229, 1843)

#plot
plot_figures(results_mat_old, programme_name = "Maternal and Elderly adult vaccination")

## cost effectiveness
table_3 <- data.frame(
  Intervention = c("Palivizimab", "Niservimab", "Maternal Vaccine","Elderly Vaccine", 
                   "Palivizimab", "Niservimab", "Maternal Vaccine","Elderly Vaccine", 
                   "Palivizimab", "Niservimab", "Maternal Vaccine","Elderly Vaccine"),
  'Cases averted' = c(16206,33064,12194,177,3159,7047,46976,177,77,319,394,477),
  'Symptomatic/GP ' = c(97,198,73,1,51,113,752,3,1,5,6,8),
  'Hosp admission ' = c(65,131,49,1,1,3,22,0,0,0,0,0), 
  'Death' = c(1,3,1,0,0,0,3,0,0,0,0,0),
  'Incremental cost' = c(" ", "", "", " ",""," ", "", "", " ","","",""),
  'Age group' = c('0-4 years','0-4 years','0-4 years','0-4 years','5-64 years','5-64 years','5-64 years','5-64 years','65+','65+','65+','65+'),
  'Prop cases (%)' = c(13.46,27.50,10.14,0.15,2.63,5.90,39.06,0.15,0.06,0.27,0.33,0.4),
  'Prop hosp (%)' = c(23.90,48.16,18.01,0.37,0.37,1.10,8.09,0,0,0,0,0),
  'Prop death'= c(12.5, 37.5,12.5,0,0,0,37.5,0,0,0,0,0),
  'Life years gained' = c(85,255,85,0,0,0,'NA',0,0,0,0,0)

)






