library(dplyr)
library(data.table)
library(readr)
library(purrr)

# Set the working directory
setwd('/Volumes/igmm/GenScotDepression/data/abcd/release5.0/core')

csv_files <- c('abcd-general/abcd_p_demo.csv', 'mental-health/mh_p_ksads_bg.csv', 
               'gender-identity-sexual-health/gish_p_gi.csv', 
               'mental-health/mh_p_fhx.csv',
               'physical-health/ph_p_sds.csv','physical-health/ph_y_pds.csv')

# Read each selected CSV file into a list of data frames using purrr
dfs <- map(csv_files, ~ read_csv(.x))

################

for (i in seq_along(dfs)) {
  df <- dfs[[i]]  # Get the current dataframe from the list
  df$eventname <- ifelse(df$eventname == 'baseline_year_1_arm_1', 0,
                         ifelse(df$eventname == '6_month_follow_up_arm_1', 1,
                                ifelse(df$eventname == '1_year_follow_up_y_arm_1', 2,
                                       ifelse(df$eventname == '18_month_follow_up_arm_1', 3,
                                              ifelse(df$eventname == '2_year_follow_up_y_arm_1', 4,
                                                     ifelse(df$eventname == '30_month_follow_up_arm_1', 5,
                                                            ifelse(df$eventname == '3_year_follow_up_y_arm_1', 6,
                                                                   ifelse(df$eventname == '42_month_follow_up_arm_1', 7, 8))))))))
  df$eventname <- as.numeric(df$eventname)
  df$src_subject_id <- as.character(df$src_subject_id)
  dfs[[i]] <- df  # Assign the modified dataframe back to the list
}

merged_df <- reduce(dfs, full_join, by = c('src_subject_id', 'eventname'), suffix = c("", ".y"))

### need to make a key due to some predictors being time invariant


# select vars
#### problem with fam_history_q6d_depression 
#variables <- c("src_subject_id", "eventname", "kbi_p_c_bully", "demo_sex_v2",
#               "fam_history_q6d_depression", "sds_p_ss_total", "demo_comb_income_v2",'pds_y_ss_female_category_2',
#               'pds_y_ss_male_category_2')

#predictors <- merged_df[variables]


################### puberty data ################
# 1 - prepuberty; 2 - early puberty; 3 - mid puberty; 4 - late puberty; 5 - post puberty

sex_puberty_vars <- c("src_subject_id", "eventname", "demo_sex_v2",
                      "pds_y_ss_female_category_2","pds_y_ss_male_cat_2")

predictors <- merged_df[sex_puberty_vars]
puberty.vars <- predictors

# fill sex from t0 downwards
puberty.vars <- puberty.vars %>%
  group_by(src_subject_id) %>%
  fill(demo_sex_v2, .direction = "down")

# Create a new column that combines the values from both columns, regardless of sex
puberty.vars$pds_y_ss <- ifelse(!is.na(puberty.vars$pds_y_ss_female_category_2), puberty.vars$pds_y_ss_female_category_2, 
                      ifelse(!is.na(puberty.vars$pds_y_ss_male_cat_2), puberty.vars$pds_y_ss_male_cat_2, NA))

# Remove the individual sex-specific columns
puberty.vars <- subset(puberty.vars, select = -c(pds_y_ss_female_category_2, pds_y_ss_male_cat_2))

write.table(puberty.vars, '/Volumes/igmm/GenScotDepression/users/poppy/abcd/predictors/puberty.vars.txt', col.names=TRUE)

################################################ 


# restrict to valid scores
predictors_qcd <- predictors %>%
  mutate(kbi_p_c_bully = ifelse(kbi_p_c_bully >= 1 & kbi_p_c_bully <= 2, kbi_p_c_bully, NA),
         demo_sex_v2 = ifelse(demo_sex_v2 >= 1 & demo_sex_v2 <= 2, demo_sex_v2, NA),
         fam_history_q6d_depression = ifelse(fam_history_q6d_depression >=0 & fam_history_q6d_depression <= 1, fam_history_q6d_depression, NA),
         sds_p_ss_total = ifelse(sds_p_ss_total >=0 & sds_p_ss_total <= 150, sds_p_ss_total, NA),
         demo_comb_income_v2 = ifelse(demo_comb_income_v2 >=1 & demo_comb_income_v2 <= 10, demo_comb_income_v2, NA),
  )

# 1 male 2 female --> recode to: 0 male 1 female
  
predictors_recoded <- predictors_qcd %>%
  mutate(kbi_p_c_bully = ifelse(kbi_p_c_bully == 1, 0, ifelse(kbi_p_c_bully == 2, 1, kbi_p_c_bully)),
         demo_sex_v2 = ifelse(demo_sex_v2 == 1, 0, ifelse(demo_sex_v2 == 2, 1, demo_sex_v2))
  )
  
  
write.table(predictors_recoded, '/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_envi_vars.txt')