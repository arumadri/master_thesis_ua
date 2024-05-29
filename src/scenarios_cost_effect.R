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
params_mat$uptake1      <- 0.6*0.46   # 0.46 = proportion of newborns in 0-4 year age group pop. (679995/ 1477479)                           
params_mat$uptake2      <- 0.6*0.02  # 0.02 = proportion of pregnant women in 5-64 year age group pop. (679995/41200378)                          
params_mat$uptake3      <- 0                                
params_mat$vac_eff1         <- 0.513                                 
params_mat$vac_eff2         <- 0.717                                   
params_mat$vac_eff3         <- 0

results_mat <- ode(y = state_initial_age, times = time_ode, func = rsv_model_age,
                    parms = params_mat, method = "ode23")

# no. of cases 
round(new_cases(results_mat), digits = 0) # 6140276 1989860   42169
## 6111777 1993124   42140 
# averted c(28499, -3264, 29)

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
table_3 <- data.frame(
  Intervention = c("Palivizumab", "Niservimab", "Maternal vaccine","Elderly vaccine", 
                   "Palivizumab", "Niservimab", "Maternal vaccine","Elderly vaccine", 
                   "Palivizumab", "Niservimab", "Maternal vaccine","Elderly vaccine"),
  'Symptomatic/GP ' = c(121,975,295,7,0,0,0,0,0,0,0,0),
  'Hosp admission ' = c(55,368,135,0,0,0,0,0,0,0,0,0), 
  'Death' = c(1,8,3,0,0,0,0,0,0,0,0,0),
  'Life years gained' = c(85,585,255,0,0,0,0,0,0,0,0,0),
  'Incr cost GP' = c(362075963, 50863410,17228266, 306805,362075963, 50863410,17228266, 306805,362075963, 50863410,17228266, 306805),
  'Incr cost hosp'= c(362038350,50443903,17074196,307021,362038350,50443903,17074196,307021,362038350,50443903,17074196,307021),
  'Age group' = c('0-4 years','0-4 years','0-4 years','0-4 years','5-64 years','5-64 years','5-64 years','5-64 years','65+','65+','65+','65+'),
  'No of cases' = c(16041,106052,38961,526,0,0,0,0,0,0,0,0),
  'No hosp' = c(55,368,135,0,0,0,0,0,0,0,0,0),
  'No death'= c(1, 8,3,0,0,0,0,0,0,0,0,0),
  'No GP' = c(121,975,295,7,0,0,0,0,0,0,0,0)
)
## no. cases averted
# palivizimab = c(13532, 2461,48)
# niservimab  = c(90270, 15548, 234)
# maternal    = c(32921, 5927, 113)
# elderly     = c(87, 57,382)

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

## averted per health outcome 
## GP 
# palivizimab = c(81, 39,1)
# niservimab  = c(542, 429, 4)
# maternal    = c(198, 95, 2)
# elderly     = c(0, 1,6)

## hospitalizations
# palivizimab = c(54, 1,0)
# niservimab  = c(361, 7, 0)
# maternal    = c(132, 3, 0)
# elderly     = c(0, 0,0)

## death
# palivizimab = c(1, 0,0)
# niservimab  = c(7, 1, 0)
# maternal    = c(3, 0, 0)
# elderly     = c(0, 0,0)

## uptake       ## 2490, 2855, 3220
# palivizimab = (97654.98+96784.53 + 96367.3)*0.9*0.338  = 88463
# niservimab  = (97654.98+96723.92 +96418.84)*0.9*0.781  = 204402
# maternal    = (96114.08+96618.84+97654.98)*0.6*0.513*0.5
# + (39507.77+38739.64+38440.74)*0.6*0.717*0.5 =  69790
# elderly     = (899.9644 + 883.8788+ 877.5723)*0.7*0.667 = 1243

## costs (admin + price )
# palivizimab = (57.50+4035.50) = 4093
# niservimab  = (11.00+238)     = 249
# maternal    = (9.00+238)      = 247
# elderly     = (9.00+238)      = 247

## GP/symptomatic = 36.00
## Hosp admin     = 1100.23 (<5yr), 652.29 (>5yr)

### incremental cost 
## GP/symptomatic 
 # palivizimab = (88463*4093) - ((56+29+1)*36.00) = 362075963
 # niservimab  = (204402*249)  - ((617+287+4)*36.00) = 50863410
 # maternal    = (69790*247) - ((224+48+2)*36.00)  = 17228266
 # elderly     = (1243*247) - ((0+1+5)*36.00)   = 306805

## Hosp admin 
## < 5 years
# palivizimab = (88463*4093)  - (37*1100.23) = 362038350
# niservimab  = (204402*249) - (411*1100.23) = 50443903
# maternal    = (69790*247) - (149*1100.23)  = 17074196
# elderly     = (1243*247)   - (0*1100.23)   = 307021

## > 5 years 
# palivizimab = (88463*4093) - (1*652.29) = 362078407
# niservimab  = (204402*249) - (8*652.29) = 50890880
# maternal    = (69790*247) - (1*652.29) = 17237478
# elderly     = (1243*247) - (0*652.29)   = 307021
