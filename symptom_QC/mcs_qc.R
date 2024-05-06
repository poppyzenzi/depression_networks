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

#### age 5 ####

#AGE AND SEX
cm_interview = read_dta("UKDA-5795-stata/stata/stata13/mcs3_cm_interview.dta")
cm_interview$ID = paste0(cm_interview$MCSID, sep = "_", cm_interview$CCNUM00)
cm_interview = cm_interview %>% filter(ID %in% IDunique) %>% select(ID, CHCSEX00, CHCAGE00)
cm_interview$CHCAGE00 = cm_interview$CHCAGE00/365 # bc age is in years

family_derived = read_dta("UKDA-5795-stata/stata/stata13/mcs3_family_derived.dta")
family_derived$ID = paste0(family_derived$MCSID, "_1")
family_derived = family_derived %>% select(ID, COECDUK0, CDOEDP00)

parent_interview = read_dta("UKDA-5795-stata/stata/stata13/mcs3_parent_interview.dta")
parent_interview = parent_interview %>% filter() %>% select(MCSID, CELIG00, CPMAFI00)

parent_unique = parent_interview[with(parent_interview, order(CELIG00)), ] %>% distinct(MCSID, .keep_all = TRUE)
parent_unique$ID = paste0(parent_unique$MCSID, "_1")

demogs_5yo = merge(cm_interview, family_derived, by = "ID")
demogs_5yo = merge(demogs_5yo, parent_unique, by = "ID", all.x = TRUE)

rm(cm_interview, family_derived, parent_interview, parent_unique)


#### age 7 ####

#AGE AND SEX
cm_interview = read_dta("UKDA-6411-stata/stata/stata13/mcs4_cm_interview.dta")
cm_interview$ID = paste0(cm_interview$MCSID, sep = "_", cm_interview$DCNUM00)
cm_interview = cm_interview %>% filter(ID %in% IDunique) %>% select(ID, DCCSEX00, DCAGEY00)
mean_age <- mean(cm_interview$DCAGEY00[cm_interview$DCAGEY00 > 1], na.rm = TRUE)

# read in dat, make unique IDs
sdq4 <- read_dta('UKDA-6411-stata/stata/stata13/mcs4_cm_derived.dta')
sdq4$ID = paste0(sdq4$MCSID, sep = "_", sdq4$DCNUM00)
IDunique = sdq4$ID[grep("_1", sdq4$ID)]
# SDQ derived variables
sdq4_selec = sdq4 %>% filter(ID %in% IDunique) %>% select(ID, DDEMOTION, DDCONDUCT, DDHYPER, DDPEER, DDPROSOC)
internalising = sdq4_selec %>% select(ID, DDEMOTION, DDPEER) %>%
  mutate(int_score = ifelse(DDEMOTION >= 0, DDEMOTION, NA) + ifelse(DDPEER >= 0, DDPEER, NA))
mean(internalising$int_score, na.rm = TRUE)
sd(internalising$int_score, na.rm = TRUE)

# symptom level data
sdq4_sym <- read_dta('UKDA-6411-stata/stata/stata13/mcs4_parent_cm_interview.dta')
sdq4_sym$ID = paste0(sdq4_sym$MCSID, sep = "_", sdq4_sym$DPNUM00)
IDunique = sdq4_sym$ID[grep("_1", sdq4_sym$ID)]
# emotional symptoms (1) and peer relationship problems (4) = internalising problems score 
sdq4_sym_selec = sdq4_sym %>% filter(ID %in% IDunique) %>% select(ID, DPSDHS00, DPSDMW00, DPSDUD00, DPSDNC00, DPSDFE00,
                                                                  DPSDSP00, DPSDGF00, DPSDLC00, DPSDPB00, DPSDGB00)

sdq4_totals = 
  
######## age 11 ######## 

#AGE
cm_interview = read_dta("UKDA-7464-stata/stata/stata13/mcs5_cm_interview.dta")
cm_interview$ID = paste0(cm_interview$MCSID, sep = "_", cm_interview$ECNUM00)
cm_interview = cm_interview %>% filter(ID %in% IDunique) %>% select(ID, EMCS5AGE)
mean_age <- mean(cm_interview$EMCS5AGE[cm_interview$EMCS5AGE > 1], na.rm = TRUE)

# read in dat, make unique IDs
sdq5 <- read_dta('UKDA-7464-stata/stata/stata13/mcs5_cm_derived.dta')
sdq5$ID = paste0(sdq5$MCSID, sep = "_", sdq5$ECNUM00)
IDunique = sdq5$ID[grep("_1", sdq5$ID)]

# SDQ derived variables
sdq5_selec = sdq5 %>% filter(ID %in% IDunique) %>% select(ID, EEMOTION, ECONDUCT, EHYPER, EPEER, EPROSOC)
internalising = sdq5_selec %>% select(ID, EEMOTION, EPEER) %>%
  mutate(int_score = ifelse(EEMOTION >= 0, EEMOTION, NA) + ifelse(EPEER >= 0, EPEER, NA))
mean(internalising$int_score, na.rm = TRUE)
sd(internalising$int_score, na.rm = TRUE)

# symptom level data
sdq5_sym <- read_dta('UKDA-7464-stata/stata/stata13/mcs5_parent_cm_interview.dta')
sdq5_sym$ID = paste0(sdq5_sym$MCSID, sep = "_", sdq5_sym$EPNUM00)
IDunique = sdq5_sym$ID[grep("_1", sdq5_sym$ID)]
# emotional symptoms (1) and peer relationship problems (4) = internalising problems score 
sdq5_sym_selec = sdq5_sym %>% filter(ID %in% IDunique) %>% select(ID, EPSDHS00, EPSDMW00, EPSDUD00, EPSDNC00, EPSDFE00,
                                                                  EPSDSP00, EPSDGF00, EPSDLC00, EPSDPB00, EPSDGB00)


#### age 14 ####

#AGE
cm_interview = read_dta("UKDA-8156-stata/stata/stata13/mcs6_cm_interview.dta")
cm_interview$ID = paste0(cm_interview$MCSID, sep = "_", cm_interview$FCNUM00)
cm_interview = cm_interview %>% filter(ID %in% IDunique) %>% select(ID, FCCAGE00)
mean_age <- mean(cm_interview$FCCAGE00[cm_interview$FCCAGE00 > 1], na.rm = TRUE)

# read in dat, make unique IDs
sdq6 <- read_dta('UKDA-8156-stata/stata/stata13/mcs6_cm_derived.dta')
sdq6$ID = paste0(sdq6$MCSID, sep = "_", sdq6$FCNUM00)
IDunique = sdq6$ID[grep("_1", sdq6$ID)]

# SDQ derived variables
sdq6_selec = sdq6 %>% filter(ID %in% IDunique) %>% select(ID, FEMOTION, FCONDUCT, FHYPER, FPEER, FPROSOC)
internalising = sdq6_selec %>% select(ID, FEMOTION, FPEER) %>%
  mutate(int_score = ifelse(FEMOTION >= 0, FEMOTION, NA) + ifelse(FPEER >= 0, FPEER, NA))
mean(internalising$int_score, na.rm = TRUE)
sd(internalising$int_score, na.rm = TRUE)

# symptom level data
sdq6_sym <- read_dta('UKDA-8156-stata/stata/stata13/mcs6_parent_cm_interview.dta')
sdq6_sym$ID = paste0(sdq6_sym$MCSID, sep = "_", sdq6_sym$FPNUM00)
IDunique = sdq6_sym$ID[grep("_1", sdq6_sym$ID)]
# emotional symptoms (1) and peer relationship problems (4) = internalising problems score 
sdq6_sym_selec = sdq6_sym %>% filter(ID %in% IDunique) %>% select(ID, FPSDHS00, FPSDMW00, FPSDUD00, FPSDNC00, FPSDFE00,
                                                                  FPSDSP00, FPSDGF00, FPSDLC00, FPSDPB00, FPSDGB00)

#### age 17 ####

#AGE
age17 = read_dta("UKDA-8682-stata/stata/stata13/mcs7_hhgrid.dta")
age17$ID = paste0(age17$MCSID, "_", age17$GCNUM00)
age17 = age17 %>% filter(ID %in% IDunique) %>% select(ID, GHCAGE00, GHCSEX00)
mean_age <- mean(age17$GHCAGE00[age17$GHCAGE00 > 1], na.rm = TRUE)

sdq7 <- read_dta('UKDA-8682-stata/stata/stata13/mcs7_cm_derived.dta')
sdq7$ID = paste0(sdq7$MCSID, sep = "_", sdq7$GCNUM00)
IDunique = sdq7$ID[grep("_1", sdq7$ID)]

# SDQ derived variables
sdq7_selec = sdq7 %>% filter(ID %in% IDunique) %>% select(ID, GEMOTION, GCONDUCT, GHYPER, GPEER, GPROSOC)
internalising = sdq7_selec %>% select(ID, GEMOTION, GPEER) %>%
  mutate(int_score = ifelse(GEMOTION >= 0, GEMOTION, NA) + ifelse(GPEER >= 0, GPEER, NA))
mean(internalising$int_score, na.rm = TRUE)
sd(internalising$int_score, na.rm = TRUE)

# symptom level data
sdq7_sym <- read_dta('UKDA-8682-stata/stata/stata13/mcs7_parent_cm_interview.dta')
sdq7_sym$ID = paste0(sdq7_sym$MCSID, sep = "_", sdq7_sym$GPNUM00)
IDunique = sdq7_sym$ID[grep("_1", sdq7_sym$ID)]
sdq7_sym_selec = sdq7_sym %>% filter(ID %in% IDunique) %>% select(ID, GPSDHS00, GPSDMW00, GPSDUD00, GPSDNC00, GPSDFE00,
                                                                  GPSDSP00, GPSDGF00, GPSDLC00, GPSDPB00, GPSDGB00)

#################################################

# get sex
demog4 <- read_dta('UKDA-6411-stata/stata/stata13/mcs4_cm_interview.dta')
demog4$ID = paste0(demog4$MCSID, sep = "_", demog4$DCNUM00)
IDunique = demog4$ID[grep("_1", demog4$ID)]
sexdf <- demog4 %>% select(ID, DCCSEX00)
# merge sex data
#sdq4_sym_selec <- merge(sdq4_sym_selec, sexdf, by='ID')

### merge 4 waves ####
# dervied dat
sdq_list <- list(sdq4_selec, sdq5_selec, sdq6_selec, sdq7_selec)
sdq_4wav <- sdq_list %>% reduce(full_join, by ="ID")
length(unique(sdq_4wav$ID)) # N=15,384
# symptoms dat
sdq_sym_list <- list(sdq4_sym_selec, sdq5_sym_selec, sdq6_sym_selec, sdq7_sym_selec)
sdq_sym_4wav <- sdq_sym_list %>% reduce(full_join, by ="ID")
length(unique(sdq_sym_4wav$ID)) # N=14958
# remove dups
sdq_symp_4wav_unique <- unique(sdq_sym_4wav)
# merge sex data
sdq_symp_4wav_unique <- merge(sdq_symp_4wav_unique, sexdf, by='ID', all.x = TRUE)



setwd('/Volumes/igmm/GenScotDepression/users/poppy/mcs')

write.table(sdq_symp_4wav_unique, 'sdq_symp_4waves_beforeQC.txt')
#################################################

## QC SYMPTOM SCORES
# -1 = NA, 1 = not true, 2 = somewhat true, 3 = certainly true, 4 = blank
# symptoms scoring -1 (NA)
# males 1, females 2 --> males 0, females 1
df <- sdq_symp_4wav_unique %>%
        replace(. == -1, NA) %>%
        replace(. == 4, NA) %>%
        replace(. == 1, 0) %>%
      replace(. == 2, 1) %>%
      replace(. == 3, 1) %>%
      replace(. == -9, NA) %>% # refusal
      replace(. == -8, NA) # don't know
# rename with sweep numbers s5 s6 s7
colnames(df) <- c("ID", 
                  "s4PSDHS00", "s4PSDMW00", "s4PSDUD00", "s4PSDNC00", "s4PSDFE00",
                  "s4PSDSP00", "s4PSDGF00","s4PSDLC00", "s4PSDPB00","s4PSDGB00",
                  
                  "s5PSDHS00", "s5PSDMW00", "s5PSDUD00", "s5PSDNC00", "s5PSDFE00",
                  "s5PSDSP00", "s5PSDGF00","s5PSDLC00", "s5PSDPB00","s5PSDGB00",
                  
                  "s6PSDHS00","s6PSDMW00", "s6PSDUD00", "s6PSDNC00", "s6PSDFE00",
                  "s6PSDSP00","s6PSDGF00", "s6PSDLC00", "s6PSDPB00", "s6PSDGB00",
                  
                  "s7PSDHS00", "s7PSDMW00", "s7PSDUD00","s7PSDNC00", "s7PSDFE00",
                  "s7PSDSP00", "s7PSDGF00", "s7PSDLC00", "s7PSDPB00", "s7PSDGB00",
                  
                  "sex"
                  )

# save wide
write.table(df, 'symptom_data/mcs_sdq_sym_wide.txt')

# first remove labels to prevent conflicting values error
df <- zap_labels(df)

# make long
df_long <- df %>%
  pivot_longer(
    cols = starts_with("s4") | starts_with("s5") | starts_with("s6") | starts_with("s7"),
    names_to = c("Sweep", ".value"),
    names_pattern = "s(\\d+)([A-Z]+.*)"
  )

# save long
write.table(df_long, 'symptom_data/mcs_sdq_sym_long.txt')

subset_df <- df_long[!is.na(rowSums(df_long[, 4:13])), ]
table(subset_df$sex, subset_df$Sweep)

