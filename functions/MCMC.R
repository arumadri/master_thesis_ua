# set up
set.seed <- seed_samples
selected_samples <- post[sample(nrow(post), 100), ]

# analysis >>> need to define model firstly 
results <- vector("list", length = nrow(selected_samples))
for (i in 1:nrow(selected_samples)) {
  
  results[[i]] <- rsv_model(params = selected_samples[i, ])
}

numeric_results <- unlist(results)

summary_stats <- summary(numeric_results)

hist(numeric_results, main = "Analysis Results", xlab = "Result")


#### params selection 
# Calculate the mean for each parameter
param_means <- colMeans(posterior_samples)

# Calculate the median for each parameter
param_medians <- apply(posterior_samples, 2, median)

# Install the HDInterval package if you haven't already
if (!requireNamespace("HDInterval", quietly = TRUE)) {
  install.packages("HDInterval")
}

library(HDInterval)

# Calculate the 95% HPD interval for each parameter
hpd_intervals <- apply(posterior_samples, 2, function(x) hdi(x, credMass = 0.95))

library(ggplot2)

# Example visualization for one parameter
ggplot(data = as.data.frame(posterior_samples), aes(x = parameter1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Posterior Distribution of Parameter 1", x = "Parameter 1", y = "Density")

# Creating a data frame for easy reporting
summary_table <- data.frame(
  Parameter = names(param_means),
  Mean = param_means,
  Median = param_medians,
  HPD_Lower = sapply(hpd_intervals, `[`, 1),
  HPD_Upper = sapply(hpd_intervals, `[`, 2)
)

print(summary_table)

