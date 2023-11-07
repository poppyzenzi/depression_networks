library(dplyr)
library(data.table)
library(readr)
library(purrr)

# Set the working directory
setwd('/Volumes/igmm/GenScotDepression/data/abcd/release5.0/core')

csv_files <- c('abcd-general/abcd_p_demo.csv', 'mental-health/mh_p_ksads_bg.csv', 
               'gender-identity-sexual-health/gish_p_gi.csv', 
               'mental-health/mh_p_fhx.csv',
               'physical-health/ph_p_sds.csv')

# Read each selected CSV file into a list of data frames using purrr
data_frames <- map(csv_files, ~ read_csv(.x))

################

# Iterate over data tables
for (df in dfs) {
    # Relabel time points and convert to numeric, and ID to character for merging
    df$eventname <- ifelse(df$eventname == 'baseline_year_1_arm_1', 0,
                     ifelse(df$eventname == '6_month_follow_up_arm_1', 0.5,
                     ifelse(df$eventname == '1_year_follow_up_y_arm_1', 1,
                     ifelse(df$eventname == '18_month_follow_up_arm_1', 1.5,
                     ifelse(df$eventname == '2_year_follow_up_y_arm_1', 2,
                     ifelse(df$eventname == '30_month_follow_up_arm_1', 2.5,
                     ifelse(df$eventname == '3_year_follow_up_y_arm_1', 3,
                     ifelse(df$eventname == '42_month_follow_up_arm_1', 3.5, 4))))))))
    df$eventname <- as.numeric(df$eventname)
    df$src_subject_id <- as.character(df$src_subject_id)
}


merged_df <- Reduce(function(x, y) merge(x, y, by = c('src_subject_id', 'eventname'), all = TRUE), dfs)
merged_df <- as.data.frame(merged_df)


# select vars
#### problem with fam_history_q6d_depression 
variables <- c("src_subject_id", "eventname", "kbi_p_c_bully", "demo_sex_v2",
               "fam_history_q6d_depression", "sds_p_ss_total", "demo_comb_income_v2")
