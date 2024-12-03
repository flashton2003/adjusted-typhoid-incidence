# Clear environment and set up
if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

# Load required libraries
library(rjags)
library(coda)
library(lubridate)
library(dplyr)
library(readr)

# Set seed for reproducibility
set.seed(012)

# Load data
qech_blood_culture_data <- read_csv("/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.03.qech_data_for_model.csv", 
                                    col_types = cols(date = col_date(format = "%Y-%m-%d")))

# Create time periods and summarize data by year
qech_blood_culture_data$year <- year(qech_blood_culture_data$date)

# Calculate counts by year
yearly_summary <- qech_blood_culture_data %>%
  group_by(year) %>%
  summarize(
    n_tested = sum(bc_tested),
    n_positive = sum(bc_positive)
  )

# Prepare data for JAGS
n_years <- nrow(yearly_summary)
n_BCpos <- yearly_summary$n_positive
n_tested <- yearly_summary$n_tested

# You'll need to add the population under surveillance for each year
# This is just a placeholder - replace with actual population data
population_by_year <- c(786363, 800264, 814164, 828065, 841966, 855867, 869768, 883668)

# Parameters for blood culture sensitivity
mu_sensitivity <- 0.59
tau_sensitivity <- 1/0.0006507705

# Define JAGS model
jags_model <- "
model {
  # Blood culture sensitivity prior
  p_BCpos ~ dnorm(mu_sensitivity, tau_sensitivity)T(0,1)
  
  # Model for all time periods
  for (t in 1:n_timeperiods) {    
    # Model for observed positive cases
    n_BCpos[t] ~ dpois(lambda_obs[t])
    lambda_obs[t] <- lambda_true[t] * p_BCpos
    
    # True cases as latent Poisson process
    true_cases[t] ~ dpois(lambda_true[t])
    
    # Very weakly informative prior for log incidence
    log(lambda_true[t]) <- beta0[t] + log(population[t])
    beta0[t] ~ dnorm(0, 1/100000000)
    
    # Derived quantities - calculated properly from observed data
    true_inc[t] <- exp(beta0[t]) * 100000  # Convert to per 100,000 person-years
    obs_inc[t] <- lambda_obs[t]/population[t] * 100000  # Use lambda_obs instead of n_BCpos directly
    adj_factor[t] <- true_inc[t]/obs_inc[t]
  }
}
"

# Prepare data for JAGS
jags_data <- list(
  n_timeperiods = n_years,
  n_BCpos = n_BCpos,
  mu_sensitivity = mu_sensitivity,
  tau_sensitivity = tau_sensitivity,
  population = population_by_year  # Add population data
)

# Parameters to monitor
params <- c("true_inc", "obs_inc", "adj_factor", "p_BCpos")

jags <- jags.model(textConnection(jags_model), 
                   data = jags_data,
                   n.chains = 3)

# Burn-in
update(jags, 10000)

# Draw samples
samples <- coda.samples(jags, 
                        variable.names = params,
                        n.iter = 50000)


#########################
###### DIAGNOSTICS ######
#########################

# Check convergence
gelman_stats <- gelman.diag(samples, multivariate = FALSE)
print(gelman_stats)

# Calculate summary statistics
summary_stats <- summary(samples)

# Extract and format results
results <- data.frame(
  year = unique(qech_blood_culture_data$year),
  observed_incidence = (yearly_summary$n_positive / population_by_year) * 100000,  # Now as rate per 100,000
  adjusted_incidence_median = summary_stats$statistics[grep("true_inc", rownames(summary_stats$statistics)), "Mean"],
  adjusted_incidence_lower = summary_stats$quantiles[grep("true_inc", rownames(summary_stats$quantiles)), "2.5%"],
  adjusted_incidence_upper = summary_stats$quantiles[grep("true_inc", rownames(summary_stats$quantiles)), "97.5%"]
)

write_csv(results, "/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.03.qech_adjusted_incidence.csv")
