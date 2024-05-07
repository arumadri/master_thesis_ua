# prepare data from results_base
incidence_model <- data.frame(
  Age_Group = c("0-4yr", "5-64yr", "65+yr"),
  Incidence = c(7687394, 2505202, 53501),  # Actual incidence numbers
  Population = 52771252  # Population at risk for each age group
)

# Calculate incidence per 100,000 population at risk
incidence_model <- incidence_model %>%
  mutate(
    Incidence_per_100k = (Incidence / Population) * 100000  # Adjusting calculation to be per 100k of at-risk population
  )

# plot
library(ggplot2)

expected_incidence <- ggplot(incidence_model, aes(x = Age_Group, y = Incidence_per_100k, fill = Age_Group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
  scale_fill_manual(values = c("0-4yr" = "red", "5-64yr" = "green", "65+yr" = "blue"), 
                    name = "",  
                    labels = c("0-4yr", "5-64yr", "65+yr")) +
  labs(title = "Burden by Age Group",
       x = "",
       y = "Annual incidence of RSV infection per 100,000",
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
incidence_base <- plot_incidence(results_base, "Burden by Age Group")

# combine plots and save 
estimated_inicidence <- (incidence_base / expected_incidence)

fig1 <- estimated_inicidence + plot_annotation(
  title = "Estimated annual burden of RSV",
  subtitle = "",
  theme = theme(
    plot.title = element_text(family="sans",size=18, color = "black", face = "bold"),
    plot.subtitle = element_text(family="sans",size=14, color = "black")+
      theme_classic()
  ))
ggsave("fig1.pdf", plot = fig1, width = 15, height = 12, dpi = 600, path = "figs/")
### 

 
