### ALSPAC script to get SDQ symptom level data from 2023 dta #####

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

# = age?
sex <- kz021

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

t01 <- c('kq343', 'kq335', 'kq332', 'kq327', 'kq322', 'kq342', 'kq338', 'kq333', 'kq330', 'kq325')
t02 <- c('tc4023', 'tc4015', 'tc4012', 'tc4007', 'tc4002', 'tc4022', 'tc4018', 'tc4013', 'tc4010', 'tc4005')
t03 <- c('ku703', 'ku695', 'ku692', 'ku687', 'ku682', 'ku702', 'ku698', 'ku693', 'ku690', 'ku685')
t04 <- c('kw6523', 'kw6515', 'kw6512', 'kw6507', 'kw6502', 'kw6522', 'kw6518', 'kw6513', 'kw6510', 'kw6505')
t05 <- c('ta7023', 'ta7015', 'ta7012', 'ta7007', 'ta7002', 'ta7022', 'ta7018', 'ta7013', 'ta7010', 'ta7005')

cols_to_select <- c("subject", "kz021", t01, t02)
alspac_23_qc <- alspac23[selected_columns]

### get age
### rename columns
### put into long format
### QC symptom scores 
  
' Emotional problems scale
ITEM 3: Often complains of headaches… (I get a lot of headaches…) 
ITEM 8: Many worries… (I worry a lot)
ITEM 13: Often unhappy, downhearted… (I am often unhappy….)
ITEM 16: Nervous or clingy in new situations… (I am nervous in new
situations…)
ITEM 24: Many fears, easily scared (I have many fears…)
Peer problems scale
ITEM 6: Rather solitary, tends to play alone (I am usually on my own)
ITEM 11: Has at least one good friend (I have one good friend or more)
ITEM 14: Generally liked by other children (Other people my age generallylike me)
ITEM 19: Picked on or bullied by other children… (Other children or young
people pick on me)
ITEM 23: Gets on better with adults than with other children (I get on better
with adults than with people my age)'