# scenarios
# 1 year
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
params_none$vac_eff1     <- 0                                 
params_none$vac_eff2     <- 0                                    
params_none$vac_eff3     <- 0

## scenario analysis 
# base: no intervention
results_base <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                        parms = params_none, method = "ode23")

# no. of cases 
round(new_cases(results_base), digits = 0) # 1621452  420636    9426 
# averted c(0,0,0)

# plot model performance over time horizon
model_stabilize <- plot_figures(results_base, programme_name = "Model stabilization", 
                                "Model becomes stable after 6 years")

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
params_pmab$vac_eff1     <- 0.338                                  
params_pmab$vac_eff2     <- 0                                    
params_pmab$vac_eff3     <- 0

results_pmab <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_pmab, method = "ode23")
# no. of cases 
round(new_cases(results_pmab), digits = 0)  # 1621452  420636    9426 
## 1620448  420570    9426 
# averted c(1004, 66,0)

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
params_nmab$vac_eff1     <- 0.781                                   
params_nmab$vac_eff2     <- 0                                    
params_nmab$vac_eff3     <- 0

results_nmab <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_nmab, method = "ode23")
# no. of cases 
round(new_cases(results_nmab), digits = 0) # 1621452  420636    9426 
## 1608494  419979    9434 
# averted c(12958, 657, -8)

#plot
niservimab_incidence <- plot_incidence(results_nmab, "Niservimab", "Black dotdash and dashed lines mark the upper and lower limit of baseline incidence respectively") 
ggsave("niservimab_incidence.pdf", plot = niservimab_incidence, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

# mat scenario [Maternal immunization]
params_mat <- params
params_mat$waning_rate1 <- 1/180                                 
params_mat$waning_rate2 <- 1/180                                 
params_mat$waning_rate3 <- 0                                  
params_mat$uptake1      <- 0.6*0.5   # 0.5 = proportion of newborns in 0-4 year age group pop. (679995/1496623)                           
params_mat$uptake2      <- 0.6*0.2  # 0.2 = proportion of pregnant women in 5-64 year age group pop. (679995/4120037.8)                          
params_mat$uptake3      <- 0                                
params_mat$vac_eff1     <- 0.513                                 
params_mat$vac_eff2     <- 0.717                                   
params_mat$vac_eff3     <- 0

results_mat <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_mat, method = "ode23")

# no. of cases 
round(new_cases(results_mat), digits = 0) # 1621452  420636    9426 
## 1617042  420677    9427 
# averted c(4410, -41, -1)

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
params_old$vac_eff1     <- 0                                   
params_old$vac_eff2     <- 0                                    
params_old$vac_eff3     <- 0.667

results_old <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                   parms = params_old, method = "ode23")
# no. of cases 
round(new_cases(results_old), digits = 0) # 1621452  420636    9426 
## 1621458  420638    9366  
# averted c(-6, -2,60) 

#plot
elderly_incidence <- plot_incidence(results_old, "Elderly vaccination", "Black dotdash and dashed lines mark the upper and lower limit of baseline incidence respectively") 
ggsave("elderly_incidence.pdf", plot = elderly_incidence, width = 15, height = 12, dpi = 600, path = "figs/supplementary/")

## cost effectiveness

## no. cases averted
# palivizimab = c(1004, 66,0)
# niservimab  = c(12958, 657, -8)
# maternal    = c(4410, -41, -1)
# elderly     = c(-6, -2,60)

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
# palivizimab = c(6, 1,0)
# niservimab  = c(78, 11, 0)
# maternal    = c(27, 0, 0)
# elderly     = c(0, 0,1)

## hospitalizations
# palivizimab = c(4, 0,0)
# niservimab  = c(52, 0, 0)
# maternal    = c(18, 0, 0)
# elderly     = c(0, 0,0)

## death
# palivizimab = c(0, 0,0)
# niservimab  = c(1, 0, 0)
# maternal    = c(0, 0, 0)
# elderly     = c(0, 0,0)

## uptake per season                     # use sum(results_base[2464, c("S01", "S11", "S21")]) to determine uptake for palivizumab and niservimab
# palivizimab = (92328)*0.9  = 83095
# niservimab  = (92328)*0.9  = 83095
# maternal    = (679995)*0.6 =  407997   # 1863*365
# elderly     = (723)*0.7    = 506       # sum(results_base[2464, c("S03", "S13", "S23")])

## costs (admin + price )
# palivizimab = (57.50+4035.50) = 4093
# niservimab  = (11.00+238)     = 249
# maternal    = (9.00+238)      = 247
# elderly     = (9.00+238)      = 247

## GP/symptomatic = 36.00
## Hosp admin     = 1100.23 (<5yr), 652.29 (>5yr)

### incremental cost 
 # palivizimab = (83095*4093) - (7*36.00) - (4*1100.23) = 340103182
 # niservimab  = (83095*249)  - (89*36.00) -(52*1100.23) - (2*652.29) = 20628934
 # maternal    = (407997*247) - (27*36.00) - (18*1100.23) = 100754483
 # elderly     = (506*247) - (1*36.00)   = 124946
