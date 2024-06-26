# packages 
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

library(pacman)

pacman::p_load(tidyverse, patchwork)

# extract incidence for 1 year from results_base
round(new_cases(results_base), digits = 0) # incidence 

# prepare data from results_base
incidence_model <- data.frame(
  Age_Group = c("0-4yr", "5-64yr", "65+yr"),
  Incidence = c(1621452, 420636, 9426)) # to estimate incidence per 100k use population at risk for each age group c(2717451, 4116508, 9502007) 

# factor age 
incidence_model <- incidence_model %>%
  mutate(Age_Group = factor(Age_Group, levels = c("65+yr", "5-64yr", "0-4yr")))

# plot
expected_incidence <- ggplot(incidence_model, aes(x = Age_Group, y = Incidence, fill = Age_Group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  scale_fill_manual(values = c("0-4yr" = "#f98f88", "5-64yr" = "#21c352", "65+yr" = "#629dff"), 
                    name = "",  
                    labels = c("0-4yr", "5-64yr", "65+yr"),
                    breaks = c("0-4yr", "5-64yr", "65+yr")) +
  labs(title = "Annual incidence by age group",
       x = "",
       y = "Annual incidence of RSV infection",
       tag = "E") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(family = "sans", size = 13, color = "black", angle = 0, hjust = 1),
    axis.text.y = element_blank(),
    legend.text = element_text(family = "sans", size = 13, color = "black"),
    legend.title = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.spacing = unit(0.5, "lines"), 
    strip.text.x = element_text(size = 13),
    plot.tag = element_text(face = "bold")
  ) + coord_flip()

# plot incidence over time horizon (year 7)
incidence_base <- plot_incidence_annual(results_base, "", "")

# combine plots and save 
estimated_inicidence <- (incidence_base / expected_incidence)

fig1 <- estimated_inicidence + plot_annotation(
  title = "Model estimated annual burden of RSV per age group (0-4 years, 5-64 years, 65+ years) by exposure level 0-3",
  subtitle = "",
  theme = theme(
    plot.title = element_text(family="sans",size=18, color = "black", face = "bold"),
    plot.subtitle = element_text(family="sans",size=14, color = "black")+
      theme_classic()
  ))
ggsave("fig1.pdf", plot = fig1, width = 15, height = 12, dpi = 600, path = "figs/")
### 
