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
# Load data
blood_culture_data <- read_csv("~/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.02.strataa_data_for_model.csv", col_types = cols(date = col_date(format = "%Y-%m-%dT%H:%M:%SZ")))

blood_culture_data <- blood_culture_data %>%
  filter(date > ymd("2016-12-31")) %>%
  filter(date < ymd('2020-01-01'))

# Create time periods
blood_culture_data$year <- cut(blood_culture_data$date, 
                                  breaks = "year",
                                  labels = FALSE)

# Calculate blood culture counts by quarter
bc_by_year <- blood_culture_data %>%
  group_by(year) %>%
  summarize(
    n_tested = sum(bc_tested),
    n_total = n(),
    bc_prob = n_tested / n_total
  )

# Prepare model inputs
n_years <- length(unique(blood_culture_data$year))
n_BCpos <- aggregate(bc_positive ~ year, 
                     data = blood_culture_data, 
                     FUN = sum)$bc_positive
n_tested <- bc_by_year$n_tested
n_total <- bc_by_year$n_total

alpha.bc <- n_tested
beta.bc <- n_total - n_tested

# Person-time at risk
persontime <- rep(100000, n_years)

# Parameters for blood culture sensitivity
mu_sensitivity <- 0.59
tau_sensitivity <- 1/0.0006507705

#########################
######### MODEL #########
#########################
jcode <-"
model{
  # Blood culture sensitivity prior
  p_BCpos ~ dnorm(mu_sensitivity, tau_sensitivity)T(0,1)


  # Model for all time periods
  for (t in 1:n_timeperiods) {    
    # Blood culture probability from beta distribution
    p_BC[t] ~ dbeta(alpha.bc[t], beta.bc[t])
    
    
    # Model for observed positive cases
    n_BCpos[t] ~ dpois(lambda_obs[t])
    lambda_obs[t] <- lambda_true[t] * p_BCpos * p_BC[t]
    
    # True incidence model
    log(lambda_true[t]) <- beta0[t] + log(persontime[t])
    beta0[t] ~ dnorm(-8, 1)  # centered around lower incidence values
    
    # Derived quantities
    true_inc[t] <- exp(beta0[t]) * 100000  # per 100,000 PY
    obs_inc[t] <- (n_BCpos[t] + 0.5)/(persontime[t]) * 100000
    adj_factor[t] <- true_inc[t]/obs_inc[t]
  }
}
"

#########################
####### RUN MODEL - quarters #######
#########################

inits <- function() {
  list(
    beta0 = rep(-8, n_years)
  )
}



# Prepare data for JAGS
jdat <- list(
  n_BCpos = n_BCpos,
  alpha.bc = alpha.bc,
  beta.bc = beta.bc,
  persontime = persontime,
  n_timeperiods = n_years,
  mu_sensitivity = mu_sensitivity,
  tau_sensitivity = tau_sensitivity
)

# Initialize model
jmod <- jags.model(textConnection(jcode), 
                   data = jdat, 
                   inits = inits,
                   n.chains = 3)

# Burn-in
update(jmod, 50000)

# Sample from posterior
jpost <- coda.samples(jmod, 
                      thin = 10,
                      c('true_inc',
                        'p_BC',
                        'p_BCpos',
                        'adj_factor'),
                      n.iter = 200000)

#########################
###### DIAGNOSTICS ######
#########################

# Check convergence
gelman_stats <- gelman.diag(jpost, multivariate = FALSE)
print(gelman_stats)

#plot(jpost)  # trace plots
#effectiveSize(jpost)  # effective sample size
#summary(jpost)

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