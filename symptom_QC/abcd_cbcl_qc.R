######### ABCD CBCL symptom QC ##########
##### for DSM-5 depression subscale #####
###### parent reported 13 symptoms ######
#########################################

library(haven)
library(tidyverse)
library(dplyr)

labels <- c("anhedonia","cries","harms selft","poor eating","worthless",
            "guilty", "tired", "sleeps less", "sleeps more", "suicidal", 
            "sleep problems", "lacks energy", "sad")

cbcl <- read.csv('/Volumes/igmm/GenScotDepression/data/abcd/release5.0/core/mental-health/mh_p_cbcl.csv')

symptoms <- c('cbcl_q05_p','cbcl_q14_p','cbcl_q18_p','cbcl_q24_p','cbcl_q35_p',
             'cbcl_q52_p', 'cbcl_q54_p', 'cbcl_q76_p', 'cbcl_q77_p', 'cbcl_q91_p',
             'cbcl_q100_p', 'cbcl_q102_p', 'cbcl_q103_p')

cols_to_select <- c('src_subject_id','eventname', symptoms)

# select cols, only for values 0-2 in symptoms (valid)
cbcl_qc <- cbcl %>% 
  subset(select = cols_to_select) %>%
  filter(if_all(all_of(symptoms), ~ . >=0 &if_all(all_of(symptoms), ~ . <=2))) %>%
  mutate(`eventname` = dplyr::recode(`eventname`,
                                     `baseline_year_1_arm_1`="0", # age 10 (not available for BPM)
                                     `6_month_follow_up_arm_1`="1", # 10.4
                                     `1_year_follow_up_y_arm_1`="2", # 10.9
                                     `18_month_follow_up_arm_1`="3", # 11.4
                                     `2_year_follow_up_y_arm_1`="4", # 12.0
                                     `30_month_follow_up_arm_1`="5", # 12.4
                                     `3_year_follow_up_y_arm_1`="6", # 12.9
                                     `42_month_follow_up_arm_1`="7", # 13.4
                                     `4_year_follow_up_y_arm_1`="8")) # 14.1


## recode scores to binary 1 true 0 not true
cbcl_binary <- cbcl_qc %>%

  mutate(across(3:ncol(.), ~case_when(
    . == 2 ~ 1,
    TRUE ~ .
  )))

length(unique(cbcl_binary$src_subject_id))
# 11,866

write.table(cbcl_binary, '/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data/cbcl_symptoms_binary.txt')



