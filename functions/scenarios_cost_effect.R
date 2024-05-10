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
# no. of cases 
round(new_cases(results_base), digits = 0) #  
# c(6140276,1989860, 42169) 
# averted c(0,0,0)

#plot
model_stabilize <- plot_figures(results_base, programme_name = "Model stabilization", 
                                "Model becomes stable after 3 years")
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
round(new_cases(results_pmab), digits = 0)  # 6140276 1989860   42169 
# c(6131020,1988066,42125)
# averted c(9256, 1794,44)


#plot
pal_plot <- plot_figures(results_pmab, programme_name = "Palivizimab for infants", 
             "Vaccination is introduced at the start of the 4th year")
ggsave("pal_plot.pdf", plot = pal_plot, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

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
round(new_cases(results_nmab), digits = 0) # 6140276 1989860   42169
# c(6037435,1971909,41897)
# averted c(102841, 17951, 272)

#plot
plot_figures(results_nmab, programme_name = "Niservimab for infants", "nmab")

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
round(new_cases(results_mat), digits = 0) # 6140276 1989860   42169
# c(6103018, 1986853, 42055)
# averted c(37258, 3007, 114)

#plot
plot_figures(results_mat, programme_name = "Maternal vaccine for mothers and infants", "mat")

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
round(new_cases(results_old), digits = 0) # 6140276 1989860   42169
# c(6140192, 1989802, 41837)
# averted c(84, 58,332) 

#plot
plot_figures(results_old, programme_name = "Elderly adult vaccination", "old")

# # # mat_old scenario [maternal & older adult vaccination]
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
#                     parms = params_mat_old, method = "ode23")
# # no. of cases averted
# round(new_cases(results_mat_old), digits = 0) # 6140276 1989860   42169
# # c(6102896, 1986783, 42153)
# # averted c(37380, 3077, 16)
# 
# #plot
# plot_figures(results_mat_old, programme_name = "Maternal and Elderly adult vaccination", "mat/old")

## cost effectiveness
table_3 <- data.frame(
  Intervention = c("Palivizimab", "Niservimab", "Maternal vaccine","Elderly vaccine", 
                   "Palivizimab", "Niservimab", "Maternal vaccine","Elderly vaccine", 
                   "Palivizimab", "Niservimab", "Maternal vaccine","Elderly vaccine"),
  'Symptomatic/GP ' = c(86,908,274,6,0,0,0,0,0,0,0,0),
  'Hosp admission ' = c(38,419,150,0,0,0,0,0,0,0,0,0), 
  'Death' = c(1,9,3,0,0,0,0,0,0,0,0,0),
  'Life years gained' = c(85,680,255,0,0,0,0,0,0,0,0,0),
  'Incr cost GP' = c(455928453, 64160259,21709587, 385104,455928453, 64160259,21709587, 385104,455928453, 64160259,21709587, 385104),
  'Incr cost hosp'= c(455890840,63740752,21555517,385320,455930449,94611198, 21718351, 385320,455930449,94611198, 21718351, 385320),
  'Age group' = c('0-4 years','0-4 years','0-4 years','0-4 years','5-64 years','5-64 years','5-64 years','5-64 years','65+','65+','65+','65+'),
  'No of cases' = c(11094,121064,40449,474,0,0,0,0,0,0,0,0),
  'Prop hosp' = c(6.3,69,24.7,0,0,0,0,0,0,0,0,0),
  'Prop death'= c(7.69, 69.23,23.08,0,0,0,0,0,0,0,0,0),
  'Prop GP' = c(6.75,71.27,21.51,0.47,0,0,0,0,0,0,0,0)
  # 'Cases averted' = c(16206,33064,12194,177,3159,7047,46976,177,77,319,394,477)
)
## no. cases averted
# palivizimab = c(9256, 1794,44)
# niservimab  = c(102841, 17951, 272)
# maternal    = c(37258, 3007, 114)
# elderly     = c(84, 58,332) 

## uptake       ## 2490, 2855, 3220
# palivizimab = (122789.6+121902.7 + 121491.2)*0.9*0.338  = 111393
# niservimab  = (122789.6+122052 +121928.8)*0.9*0.781  = 257803
# maternal    = (6078.945+5908.689+6909.682+516.0946+2987.735+7156.22) + 
#  (6079.589+5889.296+6771.534+522.2232+2985.072+6996.631) 
# + (6085.24+5887.625+6697.502+529.3+2997.567+6934.407) =  87933
# elderly     = (1124.06 + 1110.04+ 1106.243)*0.7*0.667 = 1560

## costs (admin + price )
# palivizimab = (57.50+4035.50) = 4093
# niservimab  = (11.00+238)     = 249
# maternal    = (9.00+238)      = 247
# elderly     = (9.00+238)      = 247

## GP/symptomatic = 36.00
## Hosp admin     = 1100.23 (<5yr), 652.29 (>5yr)

### incremental cost 
## GP/symptomatic 
 # palivizimab = (111393*4093) - ((56+29+1)*36.00) = 455928453
 # niservimab  = (257803*249)  - ((617+287+4)*36.00) = 64160259
 # maternal    = (87933*247) - ((224+48+2)*36.00)  = 21709587
 # elderly     = (1560*247) - ((0+1+5)*36.00)   = 385104

## Hosp admin 
## < 5 years
# palivizimab = (111393*4093)  - (37*1100.23) = 455890840
# niservimab  = (257803*249) - (411*1100.23) = 63740752
# maternal    = (87933*247) - (149*1100.23)  = 21555517
# elderly     = (1560*247)   - (0*1100.23)   = 385320

## > 5 years 
# palivizimab = (111393*4093) - (1*1100.23) = 455930449
# niservimab  = (380000*249) - (8*1100.23) = 94611198
# maternal    = (87933*247) - (1*1100.23) = 21718351
# elderly     = (1560*247) - (0*1100.23)   = 385320
