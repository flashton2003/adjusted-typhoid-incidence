# Clear environment and set up
if(!is.null(dev.list())) dev.off()
cat("\014") 
rm(list=ls())

# Load required libraries
library(rjags)
library(coda)
library(lubridate)
library(dplyr)

# Set seed for reproducibility
set.seed(012)

#########################
#### DATA LOADING #######
#########################
# Load data
blood_culture_data <- read_csv("~/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.02.strataa_data_for_model.csv", col_types = cols(date = col_date(format = "%Y-%m-%dT%H:%M:%SZ")))

blood_culture_data <- blood_culture_data %>%
  filter(date > ymd("2016-10-01")) %>%
  filter(date < ymd('2020-01-01'))

# Create time periods (quarters and years)
blood_culture_data <- blood_culture_data %>%
  mutate(
    quarter = cut(date, breaks = "quarter", labels = FALSE),
    year = year(date)
  )

# Calculate blood culture counts by quarter
bc_by_quarter <- blood_culture_data %>%
  group_by(quarter) %>%
  summarize(
    n_tested = sum(bc_tested),
    n_total = n(),
    bc_prob = n_tested / n_total,
    n_BCpos = sum(bc_positive)
  )

# Calculate blood culture counts by year
bc_by_year <- blood_culture_data %>%
  group_by(year) %>%
  summarize(
    n_tested = sum(bc_tested),
    n_total = n(),
    bc_prob = n_tested / n_total,
    n_BCpos = sum(bc_positive)
  )

# Model function to avoid code duplication
run_model <- function(n_periods, n_tested, n_total, n_BCpos, persontime_value, time_unit) {
  
  # JAGS model code
  jcode <-"
  model{
    # Blood culture sensitivity prior
    p_BCpos ~ dnorm(mu_sensitivity, tau_sensitivity)T(0,1)
    
    # Hyperpriors for beta distribution parameters
    alpha ~ dgamma(1, 0.1)
    beta ~ dgamma(1, 0.1)
    
    # Model for all time periods
    for (t in 1:n_periods) {    
      # Blood culture probability from beta distribution
      p_BC[t] ~ dbeta(alpha, beta)
      
      # Likelihood for observed blood culture data
      n_tested[t] ~ dbin(p_BC[t], n_total[t])
      
      # Model for observed positive cases
      n_BCpos[t] ~ dpois(lambda_obs[t])
      lambda_obs[t] <- lambda_true[t] * p_BCpos * p_BC[t]
      
      # True incidence model
      log(lambda_true[t]) <- beta0[t] + log(persontime[t])
      beta0[t] ~ dnorm(-8, 1)
      
      # Derived quantities
      true_inc[t] <- exp(beta0[t]) * 100000
      obs_inc[t] <- (n_BCpos[t] + 0.5)/(persontime[t]) * 100000
      adj_factor[t] <- true_inc[t]/obs_inc[t]
    }
  }
  "
  
  # Initialize model
  inits <- function() {
    list(
      p_BCpos = 0.59,
      beta0 = rep(-8, n_periods),
      alpha = 1,
      beta = 1
    )
  }
  
  # Prepare data for JAGS
  jdat <- list(
    n_BCpos = n_BCpos,
    n_tested = n_tested,
    n_total = n_total,
    persontime = rep(persontime_value, n_periods),
    n_periods = n_periods,
    mu_sensitivity = 0.59,
    tau_sensitivity = 1/0.0006507705
  )
  
  # Run model
  jmod <- jags.model(textConnection(jcode), 
                     data = jdat, 
                     inits = inits,
                     n.chains = 3)
  
  update(jmod, 50000)
  
  jpost <- coda.samples(jmod, 
                        thin = 10,
                        c('true_inc',
                          'p_BC',
                          'p_BCpos',
                          'adj_factor',
                          'alpha',
                          'beta'),
                        n.iter = 200000)
  
  # Diagnostics
  gelman_stats <- gelman.diag(jpost, multivariate = FALSE)
  print(paste("Gelman statistics for", time_unit, "analysis:"))
  print(gelman_stats)
  
  # Results
  sum_post <- summary(jpost)
  results_inc <- round(sum_post$quantiles[grep("true_inc", rownames(sum_post$quantiles)),c(3,1,5)],0)
  results_bc_prob <- round(sum_post$quantiles[grep("p_BC", rownames(sum_post$quantiles)),c(3,1,5)],3)
  
  # Save results
  write.csv(results_inc, paste0("estimated_incidence_by_", time_unit, ".csv"))
  write.csv(results_bc_prob, paste0("estimated_bc_probability_by_", time_unit, ".csv"))
  
  return(list(
    post = jpost,
    results_inc = results_inc,
    results_bc_prob = results_bc_prob,
    gelman = gelman_stats
  ))
}

#########################
####### RUN MODELS ######
#########################

# Run quarterly analysis
quarterly_results <- run_model(
  n_periods = nrow(bc_by_quarter),
  n_tested = bc_by_quarter$n_tested,
  n_total = bc_by_quarter$n_total,
  n_BCpos = bc_by_quarter$n_BCpos,
  persontime_value = 25000,
  time_unit = "quarter"
)

# Run yearly analysis
yearly_results <- run_model(
  n_periods = nrow(bc_by_year),
  n_tested = bc_by_year$n_tested,
  n_total = bc_by_year$n_total,
  n_BCpos = bc_by_year$n_BCpos,
  persontime_value = 100000,  # Adjusted for annual persontime
  time_unit = "year"
)

# Print summary comparison
cat("\nQuarterly Analysis Summary:\n")
print(summary(quarterly_results$post))

cat("\nYearly Analysis Summary:\n")
print(summary(yearly_results$post))
