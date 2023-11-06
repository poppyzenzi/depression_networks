library(tidyr)
library(dplyr)
library(haven)
library(data.table)

setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/depression_networks")

alspac <- read_dta("/Volumes/cmvm/scs/groups/ALSPAC/data/B3421_Whalley_04Nov2021.dta") 

alspac_with_id <- alspac %>%
  unite("id", c(cidB3421, qlet))

vars_to_select <- c('id', 'kz021','r2021','f8fp470','f8se126','AT5_n', 
                      'FJCI250', 'h470', 
                    #'b370', 'e390', 
                    'V1300')

predictors <- alspac_with_id %>%
  subset(select = vars_to_select)

# r2021 mother depression ever 
predictors_qcd <- predictors %>%
  mutate(kz021 = ifelse(kz021 >= 1 & kz021 <= 2, kz021, NA), # sex
         r2021 = ifelse(r2021 >= 1 & r2021 <= 3, r2021, NA), # mother dep ever
         f8fp470 = ifelse(f8fp470 >=1 & f8fp470 <= 2, f8fp470, NA), # bullying overt
         f8se126 = ifelse(f8se126 >=0 & f8se126 <= 30, f8se126, NA), # self worth / self esteem score
         AT5_n = ifelse(AT5_n >=0 & AT5_n <= 1, AT5_n, NA),  # child trauma age 0-5
         FJCI250 = ifelse(FJCI250 >=0 & FJCI250 <= 4, FJCI250, NA), # sleep symptom score
         h470 = ifelse(h470 >=1 & h470 <= 5, h470, NA), # family income per week
         V1300 = ifelse(V1300 >=1 & V1300 <= 10, V1300, NA) # monthly take home income 
  )
        

write.table(predictors_qcd, '/Volumes/igmm/GenScotDepression/users/poppy/alspac/alspac_envi_vars')
