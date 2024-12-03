library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)

# Function to create daily records from monthly counts
create_daily_records <- function(summary_file_path) {
  # Read the summary statistics Excel file
  monthly_data <- read_excel(summary_file_path)
  
  # Process the monthly data - using the datetime values directly
  monthly_processed <- monthly_data %>%
    mutate(
      # Convert directly to Date class
      date = as.Date(month),
      # Sum BC Adult and BC Paediatric
      total_cultures = `BC Adult` + `BC Paediatric`
    ) %>%
    select(date, total_cultures)
  
  # Create daily records
  daily_records <- lapply(1:nrow(monthly_processed), function(i) {
    month_date <- monthly_processed$date[i]
    count <- monthly_processed$total_cultures[i]
    
    # Get all dates in the month
    month_end <- floor_date(month_date %m+% months(1), unit = "month") - days(1)
    month_dates <- seq(month_date, month_end, by="days")
    
    # Distribute records across the month
    records_per_day <- ceiling(count / length(month_dates))
    
    # Create data frame for this month
    daily_df <- data.frame(
      date = rep(month_dates, each = records_per_day)[1:count],
      bc_tested = 1,
      bc_positive = 0
    )
    
    return(daily_df)
  })
  
  # Combine all daily records
  do.call(rbind, daily_records)
}

# Rest of the functions remain the same
load_and_filter_qech_cases <- function(filepath, organism, projects) {
  data <- read_excel(filepath, sheet = "working")
  
  filtered <- data %>%
    mutate(
      DSP = as.Date(DSP),
      year = year(DSP)
    ) %>%
    filter(
      Organism == organism,
      project %in% projects,
      `Profile Name` == "Blood Culture MC&S Procedure"
    ) %>%
    mutate(source = "script1")
  
  return(filtered)
}

process_blood_cultures <- function(summary_file_path, cases_file_path, organism, projects, output_file_path) {
  # Create base records from summary statistics
  base_records <- create_daily_records(summary_file_path)
  
  # Load positive cases
  positive_cases <- load_and_filter_qech_cases(cases_file_path, organism, projects)
  
  # For each positive case, update one record from that day
  for(case_date in unique(positive_cases$DSP)) {
    case_count <- sum(positive_cases$DSP == case_date)
    
    # Find matching records for this date
    matching_indices <- which(base_records$date == case_date & base_records$bc_positive == 0)
    
    # Update records (up to the number of cases)
    update_count <- min(case_count, length(matching_indices))
    if(update_count > 0) {
      base_records$bc_positive[matching_indices[1:update_count]] <- 1
    }
  }
  
  # Sort by date
  base_records <- base_records %>%
    arrange(date)
  
  # Write to CSV
  write.csv(base_records, output_file_path, row.names = FALSE)
  
  return(base_records)
}

# Example usage:
result <- process_blood_cultures(
  summary_file_path = "/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.11.26/BCNumbers_20241126.xlsx",
  cases_file_path = "/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.11.27/2024.11.27.iNTS_line_list.xlsx",
  organism = "Salmonella Typhimurium",
  projects = c('MLW Microbiology Diagnostics',
               'Chikwawa DHO Service',
               'Paeds Research Ward',
               'Dept of O&G'),
  output_file_path = "/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.03.qech_data_for_model.csv"
)
