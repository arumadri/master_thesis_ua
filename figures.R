plot_figures <- function(ode_output, programme_name, subtitle){
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse")
    library(tidyverse)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
    library(patchwork)
  }
  
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
       x ="",
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
       x = "",
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
       x = "",
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

############# plot incidence #############
plot_incidence <- function(ode_output, programme_name, subtitle){
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse")
    library(tidyverse)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
    library(patchwork)
  }
  
  ode_output <- as.data.frame(ode_output)
  programme_name <- as.character(programme_name)
  # 0
  results_base_stable <- as.data.frame(ode_output[2190:3650,])
  results_base_0 <- results_base_stable
  annual_incidence_0 <- results_base_0  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E01", "E02", "E03"))%>%
    mutate(state = factor(state, levels = c("E01", "E02", "E03")))
  
  annual_incidence_0 <- ggplot(annual_incidence_0, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 7380), data = subset(annual_incidence_0, state == "E01"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 446), data = subset(annual_incidence_0, state == "E02"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 2.38), data = subset(annual_incidence_0, state == "E03"), color = "black", linetype = "dashed") +
    labs(title = "Exposure level 0",
         x = "",
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
          axis.title.y = element_text(size = 13, colour = "black"),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 9, face = "bold"),
          legend.position = "none" ) 
  annual_incidence_0
  
  #1
  results_base_1 <- as.data.frame(results_base_stable)
  annual_incidence_1 <- results_base_1  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E11", "E12", "E13"))%>%
    mutate(state = factor(state, levels = c("E11", "E12", "E13")))
  
  annual_incidence_1 <- ggplot(annual_incidence_1, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 5500), data = subset(annual_incidence_1, state == "E11"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 1920), data = subset(annual_incidence_1, state == "E12"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 32), data = subset(annual_incidence_1, state == "E13"), color = "black", linetype = "dashed") +
    labs(title = "Exposure level 1",
         x = "",
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
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 9, face = "bold"),
          legend.position = "none" ) 
  
  annual_incidence_1
  # 2
  results_base_2 <- as.data.frame(results_base_stable)
  annual_incidence_2 <- results_base_2  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E21", "E22", "E23"))%>%
    mutate(state = factor(state, levels = c("E21", "E22", "E23")))
  
  annual_incidence_2 <- ggplot(annual_incidence_2, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 3820), data = subset(annual_incidence_2, state == "E21"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 2760), data = subset(annual_incidence_2, state == "E22"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 72), data = subset(annual_incidence_2, state == "E23"), color = "black", linetype = "dashed") +
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
  
  # 3
  results_base_3 <- as.data.frame(results_base_stable)
  annual_incidence_3 <- results_base_3  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E31", "E32", "E33"))%>%
    mutate(state = factor(state, levels = c("E31", "E32", "E33")))
  
  annual_incidence_3 <- ggplot(annual_incidence_3, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = 1), data = subset(annual_incidence_3, state == "E31"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 1), data = subset(annual_incidence_3, state == "E32"), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = 1), data = subset(annual_incidence_3, state == "E33"), color = "black", linetype = "dashed") +
    labs(title = "Exposure level 3",
         x = "",
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

########## plot incidence comparison #########
plot_incidence_comparison <- function(ode_output, programme_name){
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse")
    library(tidyverse)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
    library(patchwork)
  }
  
  ode_output <- as.data.frame(ode_output)
  programme_name <- as.character(programme_name)
  # 0
  results_stable <- as.data.frame(ode_output[2190:3650,])
  results_stable_0 <- results_stable
  compare_incidence_0 <- results_stable_0  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E01", "E02", "E03"))%>%
    mutate(state = factor(state, levels = c("E01", "E02", "E03")))
  
  compare_incidence_0 <- ggplot(compare_incidence_0, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Exposure level 0",
         x = "",
         y = "Population",
         tag = "A") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2190,3650), breaks = seq(2190,3650,365),
                       labels = c("0","3", "6", "9", "12") ) +
    theme(axis.text.x=element_text(family="sans",size=13, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  compare_incidence_0
  
  #1
  results_stable_1 <- as.data.frame(results_stable)
  compare_incidence_1 <- results_stable_1  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E11", "E12", "E13"))%>%
    mutate(state = factor(state, levels = c("E11", "E12", "E13")))
  
  compare_incidence_1 <- ggplot(compare_incidence_1, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Exposure level 1",
         x = "",
         y = "",
         tag = "B") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2190,3650), breaks = seq(2190,3650,365),
                       labels = c("0","3", "6", "9", "12") ) +
    theme(axis.text.x=element_text(family="sans",size=13, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  
  compare_incidence_1
  # 2
  results_stable_2 <- as.data.frame(results_stable)
  compare_incidence_2 <- results_stable_2  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E21", "E22", "E23"))%>%
    mutate(state = factor(state, levels = c("E21", "E22", "E23")))
  
  compare_incidence_2 <- ggplot(compare_incidence_2, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Exposure level 2",
         x = "Time (years)",
         y = "Population",
         tag = "C") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2190,3650), breaks = seq(2190,3650,365),
                       labels = c("0","3", "6", "9", "12") ) +
    theme(axis.text.x=element_text(family="sans",size=13, color = "black"),
          axis.text.y=element_text(family="sans",size=13, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          axis.title.x = element_text(hjust = 1, vjust = 1),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  
  compare_incidence_2
  
  # 3
  results_stable_3 <- as.data.frame(results_stable)
  compare_incidence_3 <- results_stable_3  %>% 
    gather(key = "state", value = "value", -time) %>% 
    filter(state %in% c("E31", "E32", "E33"))%>%
    mutate(state = factor(state, levels = c("E31", "E32", "E33")))
  
  compare_incidence_3 <- ggplot(compare_incidence_3, aes(x = time, y = value, color = state)) +
    geom_line() +
    geom_vline(xintercept = 2190, linetype = "dashed", color = "red") +
    labs(title = "Exposure level 3",
         x = "",
         y = "",
         tag = "D") +
    theme_minimal() +
    facet_wrap(~state, scales = "free_y") +
    theme_classic() +
    scale_x_continuous(limits = c(2190,3650), breaks = seq(2190,3650,365),
                       labels = c("0","3", "6", "9", "12") ) +
    theme(axis.text.x=element_text(family="sans",size=12, color = "black"),
          axis.text.y=element_text(family="sans",size=12, color = "black"),
          legend.text=element_text(family="sans",size=13, color = "black"),
          legend.title=element_text(family="sans",size=13, color = "black"),
          plot.tag=element_text(face="bold"),
          panel.spacing = unit(0.5, "lines"), 
          strip.text.x = element_text(size = 8),
          legend.position = "none" ) 
  
  compare_incidence_3
  
  compare_incidence <- (compare_incidence_0 | compare_incidence_1)/ (compare_incidence_2 | compare_incidence_3)
  
  
  compare_incidence <- compare_incidence + plot_annotation(
    title = programme_name,
    subtitle = "",
    theme = theme(
      plot.title = element_text(family="sans",size=16, color = "black", face = "bold"),
      plot.subtitle = element_text(family="sans",size=14, color = "black")+
        theme_classic()
    )
  )
  
  return(compare_incidence)
}

