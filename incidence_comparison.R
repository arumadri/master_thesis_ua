# extracted data from fig 1(d) Hodgson (2020) using PlotDigitizer Online App
david_model <- data.frame(incidence = c("1160778", "334191.5", "35642.61"),
                          age_group = c("0-4yr", "5-64yr", "65+yr")
                          )
david_model$incidence <- as.numeric(david_model$incidence)

# RSV model incidence for year 7 = when model has stabilized
relevant_data <- results_base[2555:2920, c("Z1", "Z2", "Z3")]

# daily changes for each age group
daily_changes <- apply(relevant_data, 2, diff)
daily_changes_df <- as.data.frame(daily_changes)

# total incidence
incidence_model <- colSums(daily_changes_df)

# save
incidence_model <- data.frame(incidence = c("1536010.2", "498134.2", "10608.3 "),
                          age_group = c("0-4yr", "5-64yr", "65+yr")
)
incidence_model$incidence <- as.numeric(incidence_model$incidence)

# combine incidences for plotting
incidence_model$Dataset <- "incidence"
david_model$Dataset <- "david"

combined_df <- bind_rows(incidence_model, david_model)

library(dplyr)
library(ggplot2)

# total population for each model
population_model <- sum(1536002.87, 497955.45, 10567.91)
population_david <- sum(1160778, 334191.5, 35642.61)

# update dataframe
incidence_model <- data.frame(
  Model = "Incidence Model",
  Incidence = c("1536002.87", "497955.45", "10567.91"),
  Age_Group = c("0-4yr", "5-64yr", "65+yr"),
  Population = population_model
)
david_model <- data.frame(
  Model = "David Model",
  Incidence = c("1160778", "334191.5", "35642.61"),
  Age_Group = c("0-4yr", "5-64yr", "65+yr"),
  Population = population_david
)

# combined updated dataframes
combined_data <- bind_rows(incidence_model, david_model) %>%
  mutate(
    Incidence = as.numeric(Incidence),  # Ensure Incidence is numeric
    Incidence_per_100k = (Incidence / Population) * 100000 
  )
# plot
library(ggplot2)

expected_incidence <- ggplot(combined_data, aes(x = Age_Group, y = Incidence_per_100k, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Estimated annual burden of RSV",
       x = "Age Group",
       y = "Annual incidence of RSV infection per 100,000",
       tag = "E") +
  scale_fill_manual(values = c("Incidence Model" = "black", "David Model" = "red"),
                    name = "",  
                    labels = c("Model", "Hodgson (2020)")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(family = "sans", size = 13, color = "black", angle = 0, hjust = 1),
    axis.text.y = element_text(family = "sans", size = 13, color = "black"),
    legend.text = element_text(family = "sans", size = 13, color = "black"),
    legend.title = element_text(family = "sans", size = 13, face = "bold", color = "black"),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.spacing = unit(0.5, "lines"), 
    strip.text.x = element_text(size = 13),
    plot.tag=element_text(face="bold")
  ) + coord_flip()

# plot incidence over time horizon (year 7)
incidence_base <- plot_incidence(results_base, "Base: No intervention")

# combine plots and save 
estimated_inicidence <- (incidence_base / expected_incidence)
ggsave("fig1.pdf", plot = estimated_inicidence, width = 15, height = 12, dpi = 600, path = "figs/")
### 

 
