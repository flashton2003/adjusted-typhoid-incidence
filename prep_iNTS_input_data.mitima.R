# Load required libraries
library(tidyverse)
library(readxl)
library(readr)


clinical_data <- read_excel("/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.11.21/MLW-Data report 2024-11-21.MITIMA blood cultures.with_residence.xlsx")


clinical_data <- clinical_data %>% filter(Residence == 'Ndirande')
clinical_data <- clinical_data %>% mutate(bc_positive = ifelse(organism == 'Salmonella Typhimurium', 1, 0))

clinical_data <- clinical_data %>% select(DSP, bc_positive)

write_csv(clinical_data, "/Users/flashton/Dropbox/GordonGroup/iNTS_Typhimurium_updates/data_munging/2024.12.02/2024.12.02.mitima_data_for_model.csv")
