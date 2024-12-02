# Clear environment and set up
if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

# Load required libraries
library(rjags)
library(coda)
library(lubridate)

# Set seed for reproducibility
set.seed(012)

#########################
#### DATA LOADING #######
#########################
# Replace this section with your real data loading
# Example structure of required data:
# blood_culture_data should be a data frame with columns:
# - date: dates of tests
# - bc_positive: binary (0/1) for typhoid positive results
# - bc_tested: binary (0/1) for whether blood culture was taken

# Example of loading real data (modify as needed):
blood_culture_data <- read_csv("~/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.02.strataa_data_for_model.csv", col_types = cols(date = col_date(format = "%Y-%m-%dT%H:%M:%SZ")))

blood_culture_data <- blood_culture_data %>%
  filter(date > ymd("2016-10-01")) %>%
  filter(date < ymd('2020-01-01'))

# Create time periods (e.g., quarters)
blood_culture_data$quarter <- cut(blood_culture_data$date, 
                                  breaks = "quarter",
                                  labels = FALSE)

# Calculate blood culture testing counts by quarter
bc_by_quarter <- aggregate(bc_tested ~ quarter, 
                           data = blood_culture_data, 
                           FUN = function(x) c(sum(x), length(x) - sum(x)))

# Prepare model inputs
n_quarters <- length(unique(blood_culture_data$quarter))
n_BCpos <- aggregate(bc_positive ~ quarter, 
                     data = blood_culture_data, 
                     FUN = sum)$bc_positive

# Person-time at risk (modify based on your population data)
persontime <- rep(250000, n_quarters) # Example: 100,000 person-years per quarter

# Parameters for blood culture sensitivity (based on literature)
mu_sensitivity <- 0.59  # mean sensitivity
tau_sensitivity <- 1/0.0006507705  # precision


#########################
######### MODEL #########
#########################
jcode <-"
model{
  # Blood culture sensitivity prior
  p_BCpos ~ dnorm(mu_sensitivity, tau_sensitivity)

  # First time period
  logit_p_BC[1] ~ dnorm(0, 0.001)  # vague prior for first period
  p_BC[1] <- 1/(1 + exp(-logit_p_BC[1]))

  # Time-varying blood culture probability for subsequent periods
  for (t in 2:n_quarters) {
    # Random walk for logit blood culture probability
    logit_p_BC[t] ~ dnorm(logit_p_BC[t-1], tau_bc)
    p_BC[t] <- 1/(1 + exp(-logit_p_BC[t]))
  }
  
  # Model for all time periods
  for (t in 1:n_quarters) {    
    # Model for observed positive cases
    n_BCpos[t] ~ dpois(lambda_obs[t])
    lambda_obs[t] <- lambda_true[t] * p_BCpos * p_BC[t]
    
    # True incidence model
    log(lambda_true[t]) <- beta0[t] + log(persontime[t])
    beta0[t] ~ dnorm(0, 1/100000000)  # weakly informative prior
    
    # Derived quantities
    true_inc[t] <- exp(beta0[t]) * 100000  # per 100,000 PY
    
    # Modified adjustment factor calculation to avoid division by zero
    obs_inc[t] <- (n_BCpos[t] + 0.5)/(persontime[t]) * 100000  # add small constant to avoid division by zero
    adj_factor[t] <- true_inc[t]/obs_inc[t]
  }
  
  # Prior for precision of random walk
  tau_bc ~ dgamma(0.001, 0.001)
}
"

#########################
####### RUN MODEL #######
#########################

# Prepare data for JAGS
jdat <- list(
  n_BCpos = n_BCpos,
  persontime = persontime,
  n_quarters = n_quarters,
  mu_sensitivity = mu_sensitivity,
  tau_sensitivity = tau_sensitivity
)

# Initialize model
jmod <- jags.model(textConnection(jcode), data=jdat, n.chains=3)

# Burn-in
update(jmod, 10000)

# Sample from posterior
jpost <- coda.samples(jmod, 
                      thin=3,
                      c('true_inc',
                        'p_BC',
                        'p_BCpos',
                        'adj_factor'),
                      n.iter=100000)

#########################
###### DIAGNOSTICS ######
#########################

# Check convergence
gelman_stats <- gelman.diag(jpost, multivariate = FALSE)
print(gelman_stats)

#########################
###### RESULTS ##########
#########################

# Extract posterior summaries
sum_post <- summary(jpost)

# Print results
results_inc <- round(sum_post$quantiles[grep("true_inc", rownames(sum_post$quantiles)),c(3,1,5)],0)
results_bc_prob <- round(sum_post$quantiles[grep("p_BC", rownames(sum_post$quantiles)),c(3,1,5)],3)

# Save results
write.csv(results_inc, "estimated_incidence_by_quarter.csv")
write.csv(results_bc_prob, "estimated_bc_probability_by_quarter.csv")