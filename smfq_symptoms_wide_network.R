library(haven)
library(tidyr)
library(dplyr)
library(purrr)
library(data.table)
library(bootnet)
library(psychonetrics)

setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/depression_networks")

# read in alspac data and make subject variable
alspac <- read_dta("/Volumes/cmvm/scs/groups/ALSPAC/data/B3421_Whalley_04Nov2021.dta") %>%
  unite("Subject", c(cidB3421, qlet))

# list of 13 smfq scores at time points 1-4
mfq_t1Vars <- c("fddp110", "fddp112", "fddp113", "fddp114", "fddp115", "fddp116", "fddp118", "fddp119", "fddp121", "fddp122", "fddp123", "fddp124", 
                "fddp125")
mfq_t2Vars <- c("ff6500", "ff6502", "ff6503", "ff6504", "ff6505", "ff6506", "ff6508", "ff6509", "ff6511", "ff6512", "ff6513", "ff6514", "ff6515")
mfq_t3Vars <- c("fg7210", "fg7212", "fg7213", "fg7214", "fg7215", "fg7216", "fg7218", "fg7219", "fg7221", "fg7222", "fg7223", "fg7224", "fg7225")
mfq_t4Vars <- c("ccs4500", "ccs4502" ,"ccs4503", "ccs4504", "ccs4505", "ccs4506", "ccs4508", "ccs4509", "ccs4511", "ccs4512", "ccs4513", "ccs4514", "ccs4515")

x <- list("mfq_t01" = mfq_t1Vars, "mfq_t02" = mfq_t2Vars, "mfq_t03" = mfq_t3Vars, "mfq_t04" = mfq_t4Vars)

vars <- c("Subject", "kz021", "c804", mfq_t1Vars, mfq_t2Vars, mfq_t3Vars, mfq_t4Vars)
# age vars "fd003c", "ff0011a", "fg0011a","ccs9991a"

# get vars needed in new frame and rename non-mfq cols
data <- alspac %>% 
  subset(select = vars) %>%
  rename("id" = "Subject", "sex" = "kz021", "ethnicity" = "c804")

data_wide <- data

col_pairs <- list(
  c("fddp110", "fddp112", "fddp113", "fddp114", "fddp115", "fddp116", "fddp118", "fddp119", "fddp121", "fddp122", "fddp123", "fddp124", "fddp125"),
  c("ff6500", "ff6502", "ff6503", "ff6504", "ff6505", "ff6506", "ff6508", "ff6509", "ff6511", "ff6512", "ff6513", "ff6514", "ff6515"),
  c("fg7210", "fg7212", "fg7213", "fg7214", "fg7215", "fg7216", "fg7218", "fg7219", "fg7221", "fg7222", "fg7223", "fg7224", "fg7225"),
  c("ccs4500", "ccs4502", "ccs4503", "ccs4504", "ccs4505", "ccs4506", "ccs4508", "ccs4509", "ccs4511", "ccs4512", "ccs4513", "ccs4514", "ccs4515")
)

new_names <- c("t1", "t2", "t3", "t4")

for (i in seq_along(col_pairs)) {
  for (j in seq_along(col_pairs[[i]])) {
    col_to_rename <- col_pairs[[i]][j]
    new_col_name <- sprintf("%s.%02d", new_names[i], j)
    data_wide <- data_wide %>%
      rename(
        !!new_col_name := !!col_to_rename
      )
  }
}

colnames(data_wide)

# this is double wide format (wide by time point and symptom)
data_wide

############# for network analysis in R ##########
x <- data_wide %>%
  gather(key, value, -id, -sex, -ethnicity) %>%
  extract(key, c("time", "question"), "(t.)\\.(..)") 

x <- as.data.frame(x)

test <- x[,c(1,4,5,6)]

# Convert data frame to a data.table
setDT(test)

# Reshape the data using dcast
reshaped <- dcast(test, id + time ~ question, value.var = "value")

print(reshaped)

#####################################

# filter for only valid scores (0-3)
smfq_symptoms <- reshaped %>%
  filter(if_all(.cols = all_of(names(reshaped)[3:15]), ~ . >= 0 & . <= 3))

# choose time point
timepoint = smfq_symptoms$time == "t1"

# run bootnet for timepoint and symptom cols 
results <- bootnet(smfq_symptoms[timepoint, 3:15], 
                   nBoots = 100, default = "ggmModSelect")

# plot network
plot(results$sample)


######### using Borsboom 2021 paper #######

# network analyses can assist in modelling temp changes
# attitudes section
# use psychonetrics 

# vars to use 
vars <- names(smfq_symptoms)[3:15]

model1 <- Ising(smfq_symptoms, vars=vars, groups="time")

####################################
############# for mplus  ###########
####################################

smfq_wide <- data_wide

# only if symptoms valid (1-3)
smfq_wide <- smfq_wide %>%
  mutate_at(vars(4:55), ~ ifelse(. %in% c(1, 2, 3, NA), ., NA))

# replace missing with -9999
smfq_wide[is.na(smfq_wide)] <- -9999

# longtiudinal measurement invariance format (LMI)
smfq_wide_LMI <- smfq_wide %>%
  mutate(across(4:55, ~case_when(
    . == 1 ~ 2,
    . == 2 ~ 1,
    . == 3 ~ 0,
    TRUE ~ .
  )))

write.table(smfq_wide_LMI, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_alspac/mplus_data/smfq_symptoms_wide.dat", col.names = FALSE)

# for RI-CLPM network (make symptoms binary)
smfq_binary <- smfq_wide %>%
  mutate(across(4:55, ~case_when(
    . == 1 ~ 1,
    . == 2 ~ 1,
    . == 3 ~ 0,
    TRUE ~ .
  )))

# checks - should be 15645 (including all missing)
length(unique(smfq_binary$id))

# drop all irrelevant cols
smfq_binary_symptoms_only <- smfq_binary[,c(4:55)]

write.table(smfq_binary_symptoms_only, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/smfq_symptoms_binary.dat", col.names = FALSE, row.names = FALSE)


