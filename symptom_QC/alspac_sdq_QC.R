### ALSPAC script to get SDQ symptom level data from 2023 dta #####
### parent reported SDQ ####
library(tidyr)
library(haven)
library(dplyr)
library(purrr)
library(data.table)

# set wd
setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/depression_networks")

# read in alspac data and make subject variable
alspac23 <- read_dta("/Volumes/cmvm/scs/groups/ALSPAC/data/B3421/B3421_Whalley_18Oct2023.dta") %>%
  unite("Subject", c(cidB3421, qlet))

## list of 10 sdq symptoms at t1-5
t01 <- c('kq322', 'kq325', 'kq327', 'kq330', 'kq332','kq333', 'kq335', 'kq338', 'kq342', 'kq343')
t02 <- c('tc4002', 'tc4005', 'tc4007', 'tc4010', 'tc4012', 'tc4013', 'tc4015', 'tc4018', 'tc4022', 'tc4023')
t03 <- c('ku682', 'ku685', 'ku687', 'ku690', 'ku692', 'ku693', 'ku695', 'ku698', 'ku702', 'ku703')
t04 <- c('kw6502', 'kw6505', 'kw6507', 'kw6510', 'kw6512', 'kw6513', 'kw6515', 'kw6518', 'kw6522', 'kw6523')
t05 <- c('ta7002', 'ta7005', 'ta7007', 'ta7010', 'ta7012', 'ta7013', 'ta7015', 'ta7018', 'ta7022', 'ta7023')

cols_to_select <- c("Subject", 'kz021', t01, t02, t03, t04, t05)
alspac_23_qc <- alspac23[cols_to_select]
col_names <- names(alspac_23_qc)

## 5 time points
col_pairs <- list(
  c('kq322', 'kq325', 'kq327', 'kq330', 'kq332','kq333', 'kq335', 'kq338', 'kq342', 'kq343'),
  c('tc4002', 'tc4005', 'tc4007', 'tc4010', 'tc4012', 'tc4013', 'tc4015', 'tc4018', 'tc4022', 'tc4023'),
  c('ku682', 'ku685', 'ku687', 'ku690', 'ku692', 'ku693', 'ku695', 'ku698', 'ku702', 'ku703'),
  c('kw6502', 'kw6505', 'kw6507', 'kw6510', 'kw6512', 'kw6513', 'kw6515', 'kw6518', 'kw6522', 'kw6523'),
  c('ta7002', 'ta7005', 'ta7007', 'ta7010', 'ta7012', 'ta7013', 'ta7015', 'ta7018', 'ta7022', 'ta7023')
)

new_names <- c("t1", "t2", "t3", "t4", "t5")

for (i in seq_along(col_pairs)) {
  for (j in seq_along(col_pairs[[i]])) {
    col_to_rename <- col_pairs[[i]][j]
    new_col_name <- sprintf("%s.%02d", new_names[i], j)
    alspac_23_qc <- alspac_23_qc %>%
      rename(
        !!new_col_name := !!col_to_rename
      )
  }
}

colnames(alspac_23_qc)

x <- alspac_23_qc %>%
  gather(key, value, -Subject, -kz021) %>%
  extract(key, c("time", "question"), "(t.)\\.(..)") %>%
  as.data.frame(.)

alspac_sdq_wide <- dcast(x, Subject + kz021 + time ~ question, value.var = "value") %>% as.data.frame(.)
print(alspac_sdq_wide)

########### cleaning and QC #############

labels <- c("malaise", 
            "worries", 
            "unhappy", 
            "anxiety", 
            "fears", 
            "solitary",
            "friends*",
            "liked*",
            "bullied",
            "adult-oriented")


alspac_sdq_wide[4:13] <- alspac_sdq_wide[4:13] %>%
  replace(. == -10, NA) %>%
  replace(. == -8, NA) %>%
  replace(. == -6, NA) %>%
  replace(. == -1, NA) %>%
  replace(. == 1, -1) %>%
  replace(. == 2, 1) %>%
  replace(. == 3, 1) %>%
  replace(. == 9, NA) 

colnames(alspac_sdq_wide) <- c('id', 'sex', 'time', unlist(labels))

# make time col numeric
mapping <- c("t1" = 1, "t2" = 2, "t3" = 3, "t4" = 4, "t5" = 5)
alspac_sdq_wide$time <- mapping[alspac_sdq_wide$time]

# save
setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
write.table(alspac_sdq_wide, file='alspac_sdq_wide_qcd.txt')

#####################################################################


'
emotional: items 3,8,13,16,24
peer: items 6,11,14,19,23
ITEM 3: Often complains of headaches… (I get a lot of headaches…) kq322
ITEM 6: Rather solitary, tends to play alone (I am usually on my own)
ITEM 8: Many worries… (I worry a lot) kq227
ITEM 11: Has at least one good friend (I have one good friend or more)
ITEM 13: Often unhappy, downhearted… (I am often unhappy….) kq332
ITEM 14: Generally liked by other children (Other people my age generallylike me)
ITEM 16: Nervous or clingy in new situations… (I am nervous in new
situations…) kq335
ITEM 19: Picked on or bullied by other children… (Other children or young
people pick on me)
ITEM 23: Gets on better with adults than with other children (I get on better
with adults than with people my age)
ITEM 24: Many fears, easily scared (I have many fears…) kq343
Peer problems scale'


# emotional problems
sdq24 <- c('kq343', 'tc4023', 'ku703', 'kw6523', 'ta7023')
sdq16 <- c('kq335', 'tc4015', 'ku695', 'kw6515', 'ta7015')
sdq13 <- c('kq332', 'tc4012', 'ku692', 'kw6512', 'ta7012')
sdq08 <- c('kq327', 'tc4007', 'ku687', 'kw6507', 'ta7007')
sdq03 <- c('kq322', 'tc4002', 'ku682', 'kw6502', 'ta7002')
# peer problems
sdq23 <- c('kq342', 'tc4022', 'ku702', 'kw6522', 'ta7022') # C D F G I 23
sdq19 <- c('kq338', 'tc4018', 'ku698', 'kw6518', 'ta7018')
sdq14 <- c('kq333', 'tc4013', 'ku693', 'kw6513', 'ta7013')
sdq11 <- c('kq330', 'tc4010', 'ku690', 'kw6510', 'ta7010')
sdq06 <- c('kq325', 'tc4005', 'ku685', 'kw6505', 'ta7005')