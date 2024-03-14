############################################## 
####### ALSPAC quality control script ########
##### get symptoms scores in wide format #####
##### for downstream network analysis ########
##############################################

library(tidyr)
library(haven)
library(dplyr)
library(purrr)
library(data.table)

# set wd
setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/depression_networks")

# read in alspac data and make subject variable
alspac <- read_dta("/Volumes/cmvm/scs/groups/ALSPAC/data/B3421_Whalley_04Nov2021.dta") %>%
  unite("Subject", c(cidB3421, qlet))

alspac23 <- read_dta("/Volumes/cmvm/scs/groups/ALSPAC/data/B3421/B3421_Whalley_18Oct2023.dta")
alspac23 <- alspac23 %>% unite("Subject", c(cidB3421, qlet))

# list of 13 smfq scores at time points 1-4
mfq_t1Vars <- c("fddp110", "fddp112", "fddp113", "fddp114", "fddp115", "fddp116", "fddp118", "fddp119", "fddp121", "fddp122", "fddp123", "fddp124", 
                "fddp125")
mfq_t2Vars <- c("ff6500", "ff6502", "ff6503", "ff6504", "ff6505", "ff6506", "ff6508", "ff6509", "ff6511", "ff6512", "ff6513", "ff6514", "ff6515")
mfq_t3Vars <- c("fg7210", "fg7212", "fg7213", "fg7214", "fg7215", "fg7216", "fg7218", "fg7219", "fg7221", "fg7222", "fg7223", "fg7224", "fg7225")
mfq_t4Vars <- c("ccs4500", "ccs4502" ,"ccs4503", "ccs4504", "ccs4505", "ccs4506", "ccs4508", "ccs4509", "ccs4511", "ccs4512", "ccs4513", "ccs4514", "ccs4515")
mfq_t5Vars <- c('CCXD900','CCXD902','CCXD903','CCXD904','CCXD905','CCXD906','CCXD908','CCXD909','CCXD911','CCXD912','CCXD913','CCXD914','CCXD915')
mfq_t6Vars <- c('cct2700','cct2701','cct2702','cct2703','cct2704','cct2705','cct2706','cct2707','cct2708','cct2709','cct2710','cct2711','cct2712')
mfq_t7Vars <- c('YPA2000', 'YPA2010', 'YPA2020', 'YPA2030', 'YPA2040', 'YPA2050', 'YPA2060', 'YPA2070','YPA2080', 'YPA2090', 'YPA2100', 'YPA2110', 'YPA2120')

x <- list("mfq_t01" = mfq_t1Vars, "mfq_t02" = mfq_t2Vars, "mfq_t03" = mfq_t3Vars, 
          "mfq_t04" = mfq_t4Vars, "mfq_t05" = mfq_t5Vars, "mfq_t06" = mfq_t6Vars, "mfq_t07" = mfq_t7Vars)

vars <- c("Subject", "kz021", "c804", mfq_t1Vars, mfq_t2Vars, mfq_t3Vars, 
          mfq_t4Vars, mfq_t5Vars, mfq_t6Vars, mfq_t7Vars)

# age vars "fd003c", "ff0011a", "fg0011a","ccs9991a"

# get vars needed in new frame and rename non-mfq cols
data <- alspac23 %>% 
  subset(select = vars) %>%
  rename("id" = "Subject", "sex" = "kz021", "ethnicity" = "c804")

data_wide <- data

col_pairs <- list(
  c("fddp110", "fddp112", "fddp113", "fddp114", "fddp115", "fddp116", "fddp118", "fddp119", "fddp121", "fddp122", "fddp123", "fddp124", "fddp125"),
  c("ff6500", "ff6502", "ff6503", "ff6504", "ff6505", "ff6506", "ff6508", "ff6509", "ff6511", "ff6512", "ff6513", "ff6514", "ff6515"),
  c("fg7210", "fg7212", "fg7213", "fg7214", "fg7215", "fg7216", "fg7218", "fg7219", "fg7221", "fg7222", "fg7223", "fg7224", "fg7225"),
  c("ccs4500", "ccs4502", "ccs4503", "ccs4504", "ccs4505", "ccs4506", "ccs4508", "ccs4509", "ccs4511", "ccs4512", "ccs4513", "ccs4514", "ccs4515"),
  c('CCXD900','CCXD902','CCXD903','CCXD904','CCXD905','CCXD906','CCXD908','CCXD909','CCXD911','CCXD912','CCXD913','CCXD914','CCXD915'),
  c('cct2700','cct2701','cct2702','cct2703','cct2704','cct2705','cct2706','cct2707','cct2708','cct2709','cct2710','cct2711','cct2712'),
  c('YPA2000', 'YPA2010', 'YPA2020', 'YPA2030', 'YPA2040', 'YPA2050', 'YPA2060', 'YPA2070','YPA2080', 'YPA2090', 'YPA2100', 'YPA2110', 'YPA2120')
  
)

new_names <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7")

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

x <- data_wide %>%
  gather(key, value, -id, -sex, -ethnicity) %>%
  extract(key, c("time", "question"), "(t.)\\.(..)") %>%
  as.data.frame(.)

selection <- x[,c(1,4,5,6)]

# Convert data frame to a data.table
setDT(selection)

# Reshape the data using dcast
reshaped <- dcast(selection, id + time ~ question, value.var = "value")

print(reshaped)
length(unique(reshaped$id))

write.table(reshaped, "/Volumes/igmm/GenScotDepression/users/poppy/alspac/smfq_symptoms_wide.txt")
write.table(reshaped, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/projects/gmm/gmm_alspac/symptom_data/smfq_symptoms_wide")

'I felt miserable or unhappy 
I didn’t enjoy anything at all 
I felt so tired I just sat around and did nothing 
I was very restless
 I felt I was no good any more 
 I cried a lot 
I found it hard to think properly and concentrate 
 I hated myself 
 I was a bad person 
I felt lonely
I thought nobody really loved me 
I thought I could never be as good as other kids 
I did everything wrong'
