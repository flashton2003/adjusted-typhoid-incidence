# Load required libraries
library(tidyverse)
library(readxl)
library(readr)

# Read in the first file
clinical_data <- read_excel("/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.11.27/ps_enrolment.cln.xlsx")

clinical_data <- clinical_data %>% select(pid, ps07_lims, ps12_dov)

# Read in the second file (lab results)
lab_results <- read_excel("/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.10.11b/MLW-Data report 2024-10-11.STRATAA isolates.xlsx") %>%
  # Filter for Salmonella Typhimurium
  filter(str_detect(tolower(organism), "salmonella typhimurium"))

lab_results <- lab_results %>% select('Request Number', 'organism')

# Perform left join
combined_data <- clinical_data %>%
  left_join(lab_results, 
            by = c("ps07_lims" = "Request Number"))

combined_data <- combined_data %>% mutate(bc_tested = ifelse(pid == 'NotRecruited', 0, 1))
combined_data <- combined_data %>% mutate(bc_positive = ifelse(organism == 'Salmonella Typhimurium', 1, 0))
combined_data$bc_positive[is.na(combined_data$bc_positive)] <- 0

combined_data <- combined_data %>% select(ps12_dov, bc_tested, bc_positive) %>% rename(date = ps12_dov)


# Optional: Save the combined dataset
write_csv(combined_data, "/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.02.strataa_data_for_model.csv")

# Preview the combined data
#glimpse(combined_data)