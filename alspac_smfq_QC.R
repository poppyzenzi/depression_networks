library(haven)
library(tidyr)
library(dplyr)
library(purrr)
library(data.table)
library(bootnet)
library(qgraph)
library(graphicalVAR)

'script to get smfq symptom scores in wide format'
  
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

setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm")

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

###############
### see other script 
# filter for only valid scores (0-3) check this? actually 1-3
smfq_symptoms <- reshaped %>%
  filter(if_all(.cols = all_of(names(reshaped)[3:15]), ~ . >= 0 & . <= 3))


smfq_binary <- smfq_symptoms %>%
  mutate(across(3:15, ~case_when(
    . == 1 ~ 1,
    . == 2 ~ 1,
    . == 3 ~ 0,
    TRUE ~ .
  )))

# checks - should be n=8787
length(unique(smfq_binary$id))








####

# estimating simplest network PMRF

pmrf <- estimateNetwork(smfq_symptoms[,3:15], 
                        default = "pcor",     #we'll check a few later
                        corMethod = "cor"     #cor, cov, spearman... #,fun
                        #, labels = c()         #a character vector if names other than columns are to be used
                        #, .dots = list()
                        #, weighted =           #logical, should the network be weighted?
                        #, signed =             #logical, should the networks' relations have signs?
                        #, directed =           #logical, is the network directed? detected automatically
)                     #Question, what are the values of the arguments weighted, signed, and directed in our case?

plot(pmrf)

#############

# list not working 
vars = as.list(colnames(smfq_symptoms[,3:15]))

graphicalVAR(smfq_symptoms, vars=vars, idvar=id, dayvar=time)


# choose time point
#timepoint = smfq_symptoms$time == "t1"

# run bootnet for timepoint and symptom cols 
results <- bootnet(symptoms[,3:15], 
                   nBoots = 100, default = "ggmModSelect")

# plot network
plot(results$sample)
