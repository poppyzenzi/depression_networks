######### symptom skew histograms ###########


bpm <- read.csv('/Volumes/igmm/GenScotDepression/data/abcd/release5.0/core/mental-health/mh_y_bpm.csv')
symptoms = c('bpm_9_y','bpm_11_y','bpm_12_y','bpm_13_y','bpm_18_y','bpm_19_y')
int_cols = c('src_subject_id','eventname', symptoms)

# list of symptom labels
labels <- c("worthless","anxious","guilty","self-conscious","unhappy","worry")

################### symptom QC ##################

# select cols for internalising, only for values 0-2 in symptoms (valid)
bpm_int <- bpm %>% 
  subset(select = int_cols) %>%
  filter(if_all(all_of(symptoms), ~ . >=0 &if_all(all_of(symptoms), ~ . <=2))) %>%
  mutate(`eventname` = dplyr::recode(`eventname`,
                                     `baseline_year_1_arm_1`="0", # age 10 (no BPM data)
                                     `6_month_follow_up_arm_1`="1", # 10.4
                                     `1_year_follow_up_y_arm_1`="2", # 10.9
                                     `18_month_follow_up_arm_1`="3", # 11.4
                                     `2_year_follow_up_y_arm_1`="4", # 12.0
                                     `30_month_follow_up_arm_1`="5", # 12.4
                                     `3_year_follow_up_y_arm_1`="6", # 12.9
                                     `42_month_follow_up_arm_1`="7", # 13.4
                                     `4_year_follow_up_y_arm_1`="8")) # 14.1


## filter for only valid symptom scores btw 0-3
bpm_symptoms <- bpm_int %>%
  filter(if_all(.cols = all_of(names(bpm_int)[3:ncol(bpm_int)]), ~ . >= 0 & . <= 3))

data_subset <- bpm_symptoms[, 3:8]
par(mfrow=c(3,3)) 

for (i in 1:ncol(data_subset)) {
  hist(data_subset[, i], breaks = -1:2, main = colnames(data_subset)[i], xlab = "Values", cex.main=0.9)
}

main_title <- "Histograms of BPM Scores in ABCD"
mtext(main_title, outer = TRUE, cex = 0.8, font = 2, line = -1)

par(mfrow=c(1, 1))
