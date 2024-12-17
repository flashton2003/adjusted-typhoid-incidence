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
strataa_blood_culture_data <- read_csv("~/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.02.strataa_data_for_model.csv", col_types = cols(date = col_date(format = "%Y-%m-%dT%H:%M:%SZ")))

strataa_blood_culture_data <- strataa_blood_culture_data %>%
  filter(date > ymd("2016-12-31")) %>%
  filter(date < ymd('2019-11-01'))

mitima_data <- read.csv("/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.02.mitima_data_for_model.plus_zing.csv", colClasses = c(DSP = "character")) %>%
  mutate(date = as.Date(DSP, format = "%d-%b-%Y")) %>% 
  select (date, bc_positive) %>% 
  mutate(bc_tested = 1) # all these have had bc_tested

# Function to add rows for missing blood cultures
add_missing_data <- function(df, monthly_data) {
  # Convert date to month-year format for matching
  df <- df %>%
    mutate(month_year = floor_date(date, "month"))
  
  # Iterate through each month in the monthly data
  for (i in 1:nrow(monthly_data)) {
    month_str <- monthly_data$date[i]
    eligible <- monthly_data$eligible[i]
    
    # Parse month-year from `monthly_data`
    start_date <- as.Date(paste0("01-", month_str), format = "%d-%b-%y")
    end_date <- start_date + days(days_in_month(start_date) - 1)
    month_year <- floor_date(start_date, "month")
    
    # Calculate the number of enrolled for the month
    enrolled <- df %>%
      filter(month_year == floor_date(start_date, "month")) %>%
      nrow()
    #print(enrolled)
    #print(eligible)
    
    # Calculate the difference
    difference <- eligible - enrolled
    
    # Generate random dates within the month
    start_date <- as.Date(paste0("01-", month_str), format = "%d-%b-%y")
    end_date <- start_date + days(days_in_month(start_date) - 1)
    random_dates <- sample(seq(start_date, end_date, by = "day"), size = difference, replace = TRUE)
    
    # Create new rows
    new_rows <- data.frame(
      date = random_dates,
      bc_tested = 0,
      bc_positive = 0
    )
    
    # Append new rows to the main dataframe
    new_rows <- new_rows %>% mutate(month_year = floor_date(date, "month"))
    df <- rbind(df, new_rows)
    # Calculate the number of enrolled for the month
    #View(df)
    #eligible <- df %>%
    #  filter(month_year == floor_date(start_date, "month")) %>%
    #  nrow()
    #print(eligible)
  }
  
  # Sort by date
  df <- df %>% arrange(date)
  return(df)
}

# Load monthly data (replace 'path_to_monthly_data.csv' with your file path)
monthly_data <- read.csv("/Users/flashton/Dropbox/GordonGroup/TyVAC/data_munging/2024.12.06/MITIMA_eligible_enrolled_num.csv") %>%
  mutate(date = as.character(date)) # Ensure date is character for processing

# need to add the eligibles but not enrolled
# calc the differene between teh enrolled (those in mitima data) and the eligibles ()
mitima_data <- add_missing_data(mitima_data, monthly_data)



mitima_data <- mitima_data %>% select(date, bc_tested, bc_positive)

mitima_data %>% mutate(month_year = floor_date(date, "month")) %>% group_by(month_year) %>% summarise(n_eligible = n(), n_tested = sum(bc_tested), n_positive = sum(bc_positive))

combined_data <- rbind(strataa_blood_culture_data, mitima_data)

# Create time periods
combined_data$year <- cut(combined_data$date, 
                                  breaks = "year",
                                  labels = FALSE)

# Calculate blood culture counts by year
bc_by_year <- combined_data %>%
  group_by(year) %>%
  summarize(
    n_tested = sum(bc_tested),
    n_total = n(),
    bc_prob = n_tested / n_total
  )

# Prepare model inputs
n_years <- length(unique(combined_data$year))
n_BCpos <- aggregate(bc_positive ~ year, 
                     data = combined_data, 
                     FUN = sum)$bc_positive
n_tested <- bc_by_year$n_tested
n_total <- bc_by_year$n_total

alpha.bc <- n_tested
beta.bc <- n_total - n_tested

# Person-time at risk
#persontime <- rep(100000, n_years)

# population estimate from script is 99476, 100197, 100918, 103081, 103802, 104523
# then, adjust for study start/end for 2019, 2022, 2024
# 2019 - until 31st October (confirm with Meiring) (0.83)
# 2022 - from May 23rd on (0.61)
# 2024 - until Nov 21st (0.89)
#persontime <- c(100000,100000,100000,58333,100000,91666)
# 2017, 18, 19, 22, 23, 24
# multiply by 0.95 because 5% of the population are more
# this one is ndiradane only
#persontime <- c(99476 * 0.915, 100197 * 0.915, (100918 * 0.83) * 0.915, (103081 * 0.61) * 0.915, 103802 * 0.915, (104523* 0.89) * 0.915)
# this one is ndirande plus zingwangwa for the mitima years
persontime <- c(99476 * 0.915, 100197 * 0.915, (100918 * 0.83) * 0.915, (153081 * 0.61) * 0.915, 153802 * 0.915, (154523* 0.89) * 0.915)


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
    # Blood culture probability from observed counts
    p_BC[t] ~ dbeta(alpha.bc[t], beta.bc[t])
    
    # Model for observed positive cases
    n_BCpos[t] ~ dpois(lambda_obs[t])
    lambda_obs[t] <- lambda_true[t] * p_BCpos * p_BC[t]
    
    # True cases as latent Poisson process
    true_cases[t] ~ dpois(lambda_true[t])
    
    # Very weakly informative prior for log incidence
    log(lambda_true[t]) <- beta0[t] + log(persontime[t])
    beta0[t] ~ dnorm(0, 1/100000000)
    
    # Derived quantities
    true_inc[t] <- exp(beta0[t]) * 100000  # per 100,000 PY
    obs_inc[t] <- (n_BCpos[t]/persontime[t]) * 100000
    adj_factor[t] <- true_inc[t]/obs_inc[t]
  }
}
"

#########################
####### RUN MODEL - years #######
#########################

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
results_inc <- as.data.frame(results_inc)
row.names(results_inc) <- c(2017, 2018, 2019, 2022, 2023, 2024)
#results_bc_prob <- round(sum_post$quantiles[grep("p_BC", rownames(sum_post$quantiles)),c(3,1,5)],3)

# Save results
write.csv(results_inc, "/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.03.ndirande_adjusted_incidence.with_zing.csv", row.names = TRUE)
#write.csv(results_bc_prob, "estimated_bc_probability_by_quarter.csv")