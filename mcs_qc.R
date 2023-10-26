#################################################################
# Millenium Cohort Study Quality Control 
# s5 = age 11, s6 = 14, s7 = 17
#################################################################

'1-4 = total difficulties score
1 & 4 = internalising
2 & 3 = externalising'

library(haven)
library(tidyverse)
library(dplyr)

setwd('/Volumes/igmm/GenScotDepression/users/poppy/mcs')

#### age 11 ####
# read in dat, make unique IDs
sdq5 <- read_dta('UKDA-7464-stata/stata/stata13/mcs5_cm_derived.dta')
sdq5$ID = paste0(sdq5$MCSID, sep = "_", sdq5$ECNUM00)
IDunique = sdq5$ID[grep("_1", sdq5$ID)]
# SDQ derived variables
sdq5_selec = sdq5 %>% filter(ID %in% IDunique) %>% select(ID, EEMOTION, ECONDUCT, EHYPER, EPEER, EPROSOC)
# symptom level data
sdq5_sym <- read_dta('UKDA-7464-stata/stata/stata13/mcs5_parent_cm_interview.dta')
sdq5_sym$ID = paste0(sdq5_sym$MCSID, sep = "_", sdq5_sym$EPNUM00)
IDunique = sdq5_sym$ID[grep("_1", sdq5_sym$ID)]
# emotional symptoms (1) and peer relationship problems (4) = internalising problems score 
sdq5_sym_selec = sdq5_sym %>% filter(ID %in% IDunique) %>% select(ID, EPSDHS00, EPSDMW00, EPSDUD00, EPSDNC00, EPSDFE00,
                                                                  EPSDSP00, EPSDGF00, EPSDLC00, EPSDPB00, EPSDGB00)

#### age 14 ####
# read in dat, make unique IDs
sdq6 <- read_dta('UKDA-8156-stata/stata/stata13/mcs6_cm_derived.dta')
sdq6$ID = paste0(sdq6$MCSID, sep = "_", sdq6$FCNUM00)
IDunique = sdq6$ID[grep("_1", sdq6$ID)]
# SDQ derived variables
sdq6_selec = sdq6 %>% filter(ID %in% IDunique) %>% select(ID, FEMOTION, FCONDUCT, FHYPER, FPEER, FPROSOC)
# symptom level data
sdq6_sym <- read_dta('UKDA-8156-stata/stata/stata13/mcs6_parent_cm_interview.dta')
sdq6_sym$ID = paste0(sdq6_sym$MCSID, sep = "_", sdq6_sym$FPNUM00)
IDunique = sdq6_sym$ID[grep("_1", sdq6_sym$ID)]
# emotional symptoms (1) and peer relationship problems (4) = internalising problems score 
sdq6_sym_selec = sdq6_sym %>% filter(ID %in% IDunique) %>% select(ID, FPSDHS00, FPSDMW00, FPSDUD00, FPSDNC00, FPSDFE00,
                                                                  FPSDSP00, FPSDGF00, FPSDLC00, FPSDPB00, FPSDGB00)

#### age 17 ####
sdq7 <- read_dta('UKDA-8682-stata/stata/stata13/mcs7_cm_derived.dta')
sdq7$ID = paste0(sdq7$MCSID, sep = "_", sdq7$GCNUM00)
IDunique = sdq7$ID[grep("_1", sdq7$ID)]
sdq7_selec = sdq7 %>% filter(ID %in% IDunique) %>% select(ID, GEMOTION, GCONDUCT, GHYPER, GPEER, GPROSOC)
# symptom level data
sdq7_sym <- read_dta('UKDA-8682-stata/stata/stata13/mcs7_parent_cm_interview.dta')
sdq7_sym$ID = paste0(sdq7_sym$MCSID, sep = "_", sdq7_sym$GPNUM00)
IDunique = sdq7_sym$ID[grep("_1", sdq7_sym$ID)]
sdq7_sym_selec = sdq7_sym %>% filter(ID %in% IDunique) %>% select(ID, GPSDHS00, GPSDMW00, GPSDUD00, GPSDNC00, GPSDFE00,
                                                                  GPSDSP00, GPSDGF00, GPSDLC00, GPSDPB00, GPSDGB00)

#################################################

### merge 3 waves ####
# dervied dat
sdq_list <- list(sdq5_selec, sdq6_selec, sdq7_selec)
sdq_3wav <- sdq_list %>% reduce(full_join, by ="ID")
length(unique(sdq_3wav$ID)) # N=14,150
# symptoms dat
sdq_sym_list <- list(sdq5_sym_selec, sdq6_sym_selec, sdq7_sym_selec)
sdq_sym_3wav <- sdq_sym_list %>% reduce(full_join, by ="ID")
length(unique(sdq_sym_3wav$ID)) # N=13,558
# remove dups
sdq_symp_3wav_unique <- unique(sdq_sym_3wav)

#################################################


# symptoms scoring -1 (NA), 1-4 likert 
df <- sdq_symp_3wav_unique %>%
        replace(. == -1, NA)

# rename with sweep numbers s5 s6 s7
colnames(df) <- c("ID", "s5PSDHS00", "s5PSDMW00", "s5PSDUD00", "s5PSDNC00", "s5PSDFE00", "s5PSDSP00", "s5PSDGF00",
                  "s5PSDLC00", "s5PSDPB00", "s5PSDGB00", "s6PSDHS00", "s6PSDMW00", "s6PSDUD00", "s6PSDNC00", "s6PSDFE00",
                  "s6PSDSP00", "s6PSDGF00", "s6PSDLC00", "s6PSDPB00", "s6PSDGB00", "s7PSDHS00", "s7PSDMW00", "s7PSDUD00",
                  "s7PSDNC00", "s7PSDFE00", "s7PSDSP00", "s7PSDGF00", "s7PSDLC00", "s7PSDPB00", "s7PSDGB00"
                  )

# fix reverse scored items ?

# save wide
write.table(df, 'symptom_data/mcs_sdq_sym_wide.txt')

# make long
df_long <- df %>%
  pivot_longer(
    cols = starts_with("s5") | starts_with("s6") | starts_with("s7"),
    names_to = c("Sweep", ".value"),
    names_pattern = "s(\\d+)([A-Z]+.*)"
  )

# save long
write.table(df_long, 'symptom_data/mcs_sdq_sym_long.txt')
