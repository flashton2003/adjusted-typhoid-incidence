# Load required libraries
library(tidyverse)
library(readxl)
library(readr)


# Load the data
clinical_data <- read_excel("/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.11.21/MLW-Data report 2024-11-21.MITIMA blood cultures.with_residence.xlsx")

# Filter rows where Profile Name is 'Blood Culture MC&S Procedure' into bc_mcs
bc_mcs <- clinical_data %>%
  filter(`Profile Name` == 'Blood Culture MC&S Procedure')

# Keep other rows in a separate dataframe called output_data
output_data <- clinical_data %>%
  filter(`Profile Name` != 'Blood Culture MC&S Procedure') %>%
  mutate(bc_positive = 0) %>%
  select(DSP, bc_positive)

# Get dates with Salmonella Typhimurium in bc_mcs
dates_with_salmonella <- bc_mcs %>%
  filter(organism == 'Salmonella Typhimurium') %>%
  pull(DSP)

# Update bc_positive for one row per date in dates_with_salmonella
output_data <- output_data %>%
  group_by(DSP) %>%
  mutate(bc_positive = ifelse(DSP %in% dates_with_salmonella & row_number() == 1, 1, bc_positive)) %>%
  ungroup()


write_csv(output_data, '/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.02.mitima_data_for_model.plus_zing.csv')

#output_data %>% mutate() group_by() %>% summarise(n = n())
