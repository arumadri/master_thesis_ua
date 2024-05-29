# packages

if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)

pacman::p_load(tidyverse, viridisLite, viridis, patchwork)

# save current table 3 as table 3_1 to avoid losing data
### save so it doesnt get overridden by table_3
write.csv(table_3, file = "table_3_1.csv")

# prepare 
long_data_3 <- table_3 %>%
  select(Intervention, No.GP, No.hosp, No.death) %>%
  gather(key = "Metric", value = "Value", -Intervention)

# plot proportions
condition_plot <- ggplot(long_data_3, aes(x = Intervention, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_fill_manual(values = c("No.GP" = "#fde724", 
                               "No.hosp" = "#440154", 
                               "No.death" = "#20908c"),
                    labels = c("No.GP" = "GP visit", "No.hosp" = "Hospitalization", "No.death" = "Death")) +
  labs(title = "Number of health outcomes averted per intervention",
       x = "", y = "Percentage (%)", tag = "B") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(family = "sans", size = 12, color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(family = "sans", size = 12, color = "black"),
        legend.text = element_text(family = "sans", size = 12, color = "black"),
        legend.title = element_blank(),  # removes the legend title
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10),
        plot.background = element_rect(fill = "white", colour = "black"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0.5, "lines"),
        plot.tag=element_text(face="bold"),
        strip.text.x = element_text(size = 5)) 

# plot number of cases by intervention
table_3_arrange <- table_3 %>%
  arrange(desc(No.of.cases)) %>%
  mutate(Intervention = factor(Intervention, levels = unique(Intervention)))

cases_averted <- ggplot(table_3_arrange, aes(x = Intervention, y = No.of.cases, fill = Intervention)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  scale_fill_manual(values = c("Palivizumab" = "#4169E1", 
                               "Niservimab" = "#000080", 
                               "Maternal vaccine" = "#1f77b4", 
                               "Elderly vaccine" = "#00BFFF")) +
  labs(title = "Number of cases averted per intervention",
       x = "", y = "",
       tag = "A") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "sans", size = 12, color = "black"),
        legend.text = element_text(family = "sans", size = 12, color = "black"),
        legend.title = element_blank(),  
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10),
        plot.background = element_rect(fill = "white", colour = "black"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 5),
        plot.tag=element_text(face="bold")) 

# plot costs 
table_3_2 <- table_3 %>%
  # filter(Intervention != "Elderly vaccine") %>%
  mutate(
    Cost_per_GP_visit = ifelse(Symptomatic.GP. > 0, Incr.cost.GP / Symptomatic.GP., 0),
    Cost_per_Hospitalization = ifelse(Hosp.admission. > 0, Incr.cost.hosp / Hosp.admission., 0),
    Cost_per_Death = ifelse(Death > 0, 362075963 / Death, 0),  
    Cost_per_Life_Year_Gained = c(4259717,86946,67561.83,0,0,0,0,0,0,0,0,0))


# transform the costs 
cost_data <- table_3_2 %>%
  select(Intervention, Cost_per_GP_visit, Cost_per_Hospitalization, Cost_per_Death, Cost_per_Life_Year_Gained) %>%
  pivot_longer(cols = starts_with("Cost"), names_to = "Category", values_to = "Cost") %>%
  mutate(Category = recode(Category,
                           'Cost_per_GP_visit' = 'GP visit',
                           'Cost_per_Hospitalization' = 'Hospitalization',
                           'Cost_per_Death' = 'Death',
                           'Cost_per_Life_Year_Gained' = 'Life years gained'
  )) %>%
  filter(Category != "Death") 

# plot
incr_costs <- ggplot(cost_data, aes(x = Intervention, y = Cost, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  scale_fill_manual(values = c("GP visit" = "#fde724", 
                               "Hospitalization" = "#440154", 
                               "Life years gained" = "#20908c")) +
  labs(title = "Incremental costs of interventions per health outcome",
       x = "", y = "Cost (£GBP)", tag = "C") +
  geom_hline(yintercept = 20000, linetype = "dashed", color = "black") + 
  annotate("text", x = 1.38, y = 3.5e6, label = "£20,000 (willingness-to-pay threshold)", vjust = -0.5, color = "black", size = 4) +
  geom_segment(aes(x = 1.5, y = 2.5e6, xend = 1.5, yend = 20000), 
               arrow = arrow(type = "closed", length = unit(0.15, "inches")), color = "black") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_text(family = "sans", size = 12, color = "black"),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(family = "sans", size = 12, color = "black"),
        axis.text.y = element_text(family = "sans", size = 12, color = "black"),
        legend.text = element_text(family = "sans", size = 12, color = "black"),
        legend.title = element_blank(),  
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10),
        plot.background = element_rect(fill = "white", colour = "black"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 5),
        plot.tag=element_text(face="bold")) + coord_flip()

figure_3 <- (cases_averted | condition_plot) / (incr_costs)

fig_3 <- figure_3 + plot_annotation(
  title = "Impact and health economics of interventions",
  subtitle = "",
  theme = theme(
    plot.title = element_text(family="sans",size=14, color = "black", face = "bold"),
    plot.subtitle = element_text(family="sans",size=14, color = "black")+
      theme_classic()
  ))

fig_3

ggsave("fig3.pdf", plot = fig_3, width = 15, height = 12, dpi = 600, path = "figs/")
