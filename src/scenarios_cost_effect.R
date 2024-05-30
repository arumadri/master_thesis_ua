# scenarios
# try 1 year
############################ run scenarios #####################

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)

pacman::p_load(deSolve, tidyverse, knitr, kableExtra, grid, gridExtra)

# time steps
time_ode <- seq(0, 3650, 1)

# model run
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

results_base_run <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_none, method = "ode23")

# plot model performance over time horizon
model_stabilize <- plot_figures(results_base_run, programme_name = "Model stabilization", 
                                "Model becomes stable after 6 years")

## scenario analysis 
# base: no intervention
results_base <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                        parms = params_none, method = "ode23")

# no. of cases 
round(new_cases(results_base), digits = 0) #  
# c(6140276 1989860   42169 ) 
# averted c(0,0,0)

# plot incidence
baseline_incidence <- plot_incidence(results_base, "Baseline: No intervention", "Black dotdashed and dashed lines mark the upper and lower limit of baseline incidence respectively")

# save 
ggsave("fig_1.pdf", plot = model_stabilize, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")
ggsave("baseline_incidence.pdf", plot = baseline_incidence, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

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
## 6126744 1987399   42121 
# averted c(13532, 2461,48)

#plot
palivizumab_incidence <- plot_incidence(results_pmab, "Palivizumab", "Black dotdash and dashed lines mark the upper and lower limit of baseline incidence respectively") 
ggsave("palivizumab_incidence.pdf", plot = palivizumab_incidence, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

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
## 6050006 1974312   41935 
# averted c(90270, 15548, 234)

#plot
niservimab_incidence <- plot_incidence(results_nmab, "Niservimab", "Black dotdash and dashed lines mark the upper and lower limit of baseline incidence respectively") 
ggsave("niservimab_incidence.pdf", plot = niservimab_incidence, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

# mat scenario [Maternal immunization]
params_mat <- params
params_mat$waning_rate1 <- 1/180                                 
params_mat$waning_rate2 <- 1/180                                 
params_mat$waning_rate3 <- 0                                  
params_mat$uptake1      <- 0.6*0.5   # 0.5 = proportion of newborns in 0-4 year age group pop. (679995/ 1477479)                           
params_mat$uptake2      <- 0.6*0.2  # 0.2 = proportion of pregnant women in 5-64 year age group pop. (679995/41200378)                          
params_mat$uptake3      <- 0                                
params_mat$vac_eff1         <- 0.513                                 
params_mat$vac_eff2         <- 0.717                                   
params_mat$vac_eff3         <- 0

results_mat <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_mat, method = "ode23")

# no. of cases 
round(new_cases(results_mat), digits = 0) # 6140276 1989860   42169
## 6108710 1989539   42106
# averted c(31566, 321, 63)

#plot
maternal_incidence <- plot_incidence(results_mat, "Maternal vaccination", "Black dotdash and dashed lines mark the upper and lower limit of baseline incidence respectively") 
ggsave("maternal_incidence.pdf", plot = maternal_incidence, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

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
## 6140189 1989803   41787 
# averted c(87, 57,382) 

#plot
elderly_incidence <- plot_incidence(results_old, "Elderly vaccination", "Black dotdash and dashed lines mark the upper and lower limit of baseline incidence respectively") 
ggsave("elderly_incidence.pdf", plot = elderly_incidence, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

## cost effectiveness

## no. cases averted
# palivizimab = c(13532, 2461,48)
# niservimab  = c(90270, 15548, 234)
# maternal    = c(31566, 321, 63)
# elderly     = c(87, 57,382)

# no. cases averted per season
# palivizimab = c(4511, 820,16)
# niservimab  = c(30090, 5183, 78)
# maternal    = c(10522, 107, 21)
# elderly     = c(29, 19,127)

## per-infection probabilities to estimate clinical outcomes 
## averted symptomatic/GP probability per-infection 
# <5 years: 0.006
# ≥5 years: 0.016

## averted hospitalization probability per-infection 
# <5 years: 0.004
# 5-64 years: 4.688 x 10^-5
# 65+ years: 6.197 x 10^-5

## averted death probability per-infection 
# <5 years: 8.197 x 10^-6
# ≥5 years: 5.697 x 10^-6

## averted per health outcome per season
## GP 
# palivizimab = c(27, 13,0)
# niservimab  = c(181, 143, 1)
# maternal    = c(63, 1, 0)
# elderly     = c(0, 0,2)

## hospitalizations
# palivizimab = c(18, 0,0)
# niservimab  = c(120, 2, 0)
# maternal    = c(42, 1, 0)
# elderly     = c(0, 0,0)

## death
# palivizimab = c(0, 0,0)
# niservimab  = c(2, 0, 0)
# maternal    = c(1, 0, 0)
# elderly     = c(0, 0,0)

## uptake per season  
# palivizimab = (97654.98)*0.9*0.338  = 29706
# niservimab  = (97654.98)*0.9*0.781  = 68641
# maternal    = (679995)*0.6 =  407997
# elderly     = (899.9644)*0.7*0.667 = 420

## costs (admin + price )
# palivizimab = (57.50+4035.50) = 4093
# niservimab  = (11.00+238)     = 249
# maternal    = (9.00+238)      = 247
# elderly     = (9.00+238)      = 247

## GP/symptomatic = 36.00
## Hosp admin     = 1100.23 (<5yr), 652.29 (>5yr)

### incremental cost 
## GP/symptomatic 
 # palivizimab = (29706*4093) - (40*36.00) = 121585218
 # niservimab  = (68641*249)  - (325*36.00) = 17079909
 # maternal    = (407997*247) - (64*36.00)  = 100772955
 # elderly     = (420*247) - (2*36.00)   = 103668

## Hosp admin 
## < 5 years
# palivizimab = (29706*4093)  - (18*1100.23) = 121566854
# niservimab  = (68641*249) - (120*1100.23) = 16959581
# maternal    = (407997*247) - (42*1100.23)  = 100729049
# elderly     = (420*247)   - (0*1100.23)   = 103740

## > 5 years 
# palivizimab = (29706*4093) - (0*652.29) = 121586658
# niservimab  = (68641*249) - (2*652.29) = 17090304
# maternal    = (407997*247) - (0*652.29) = 100775259
# elderly     = (420*247) - (0*652.29)   = 103740