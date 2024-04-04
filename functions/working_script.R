library(HDInterval)
library(tidyverse)

# 95% HPD interval for each parameter
hpd_intervals <- apply(post, 2, function(x) hdi(x, credMass = 0.95))

# Example visualization for one parameter
ggplot(data = as.data.frame(post), aes(x = xi)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Posterior Distribution of Parameter 1", x = "Parameter 1", y = "Density")
