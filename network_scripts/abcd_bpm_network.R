library(dplyr)
library(qgraph)
library(tidyr)
library(bootnet)
library(gridExtra)
library(mice)
library(VIM)

# split into QC and network scripts 

# read in symptom data
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


# recode scores to binary 1 true 0 not true
bpm_binary <- bpm_int %>%
  mutate(across(3:ncol(.), ~case_when(
    . == 2 ~ 1,
    TRUE ~ .
  )))

# check internalising symptom df
bpm_int <- bpm_binary
head(bpm_int)

# save df (long)
setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")
write.table(bpm_int, 'abcd_bpm_sym_long.txt')


################### symptom tables ################### 

## frequency table
freq_dat <- bpm_int
colnames(freq_dat) <- c("id","time", unlist(labels))

frequency_table <- freq_dat %>%
  gather(symptom, value, `worthless`:`worry`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time/sweep
  summarize(count = sum(value)) %>%
  spread(time, count, fill = 0)

print(frequency_table)
write.csv(frequency_table, file="abcd_symptom_frequencies.csv")

## proportion table
prop_table <- freq_dat %>%
  gather(symptom, value, `worthless`:`worry`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time
  summarize(proportion = sum(value == 1) / n()) %>% # calculate proportion of "yes" responses
  spread(time, proportion, fill = 0)

print(prop_table)
write.csv(prop_table, file="abcd_symptom_proportions.csv")


################### Begin symptom network analysis ################### 

## Step 1: Estimate network

theme = 'colorblind'

# create empty list
network_list <- list()

# estimate network and loop over waves, use ggmModSelect
for (wave in c(2,4,6,8)) {
  timepoint <- bpm_int$eventname == wave
  data_subset <- bpm_int[timepoint, 3:ncol(bpm_int)]
  network <- estimateNetwork(data_subset,default="ggmModSelect")
  network_list[[as.character(wave)]] <- network
}

# Loop through the results and generate plot for estimated networks
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  qgraph(results$graph, layout='spring', theme=theme, labels = labels)
  title(paste("Wave =", wave), adj=0.8, line= -0.9)
  # optional plot centrality estimate results
  #centralityPlot(results, include = c("Betweenness","Closeness","Strength"), labels=labels) 
}

## Step 2: Run bootnet resampled network

# create empty list
boot1_result_list <- list()

# for each wave run bootnet and store results in list 
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  boot1<-bootnet(results, nBoots=100, default="ggmModSelect", nCores=1)
  boot1_result_list[[as.character(wave)]] <- boot1
}

# plot bootstrapped network results for each wave
for (wave in names(boot1_result_list)) {
  results <- boot1_result_list[[wave]]
  plot(results$sample, label = labels, theme='colorblind')
  title(paste("Wave =", wave), adj=0.8, line= -0.9) 
}

# plot bootstrapped edge CIs
# Create a function for list of plots and arrange in 1xn grid
plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, labels = FALSE, order = "sample")})
# plot
gridExtra::grid.arrange(grobs = plot_list, ncol = 4)

# Create a list of edge (bootstrapped difference tests btw non-zero edge weights) difference plots
edge_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample", labels=FALSE)})
# plot significant differences (alpha=0.05) of edges
# grey = non sig diff, black = sig diff, coloured = colour of edge
grid.arrange(grobs = edge_plot_list, ncol = 4)

# Create a list of strength difference plots and arrange in 1xn grid
strength_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "strength", plot = "difference", order = "sample")})
# Plot significant differences (alpha=0.05) of nodestrength
# grey = non sig, black = sig, white = value of node strength 
grid.arrange(grobs = strength_plot_list, ncol = 4)

# test for difference in strength between two nodes
node1 = "bpm_19_y"
node2 = "bpm_12_y"
differenceTest(boot1, node1, node2, "strength")

################### Adding environmental data ################### 

# read in environmental vars, merge with symptom data, QC var names
predictors <- read.table('/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_envi_vars.txt')

# predictors recorded at baseline (time invariant)
t0_predictors <- predictors %>% subset(predictors$eventname==0) %>% select(-c('eventname'))
bpm_env <- merge(bpm_int, t0_predictors, by=c('src_subject_id')) %>% select(-c('sds_p_ss_total'))

# longitudinal predictors (time variant), here sleep only recorded t0-4, bpm we dont have t0, so just t1-t4
bpm_env <- merge(bpm_env, predictors[c("src_subject_id", "eventname", "sds_p_ss_total")], 
                 by = c("src_subject_id", "eventname"), 
                 all.x = TRUE) 

# rename cols
bpm_env <- bpm_env %>% rename(`sex` = demo_sex_v2,
                                `mother_depression` = fam_history_q6d_depression,
                                `bullying` = kbi_p_c_bully,
                                #`child trauma` = ,
                                `sleep` = sds_p_ss_total,
                                `income` = demo_comb_income_v2,
)

# get environmental col labels
env_labels = list(names(bpm_env[,9:ncol(bpm_env)]))

################### Adding genetic data ################### 

# read in genetic vars, merge with symptom data, QC var names
setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/projects/prs/prs_bpm_OUT')
prs <- read.table('all_thresholds/abcd_mddipsych_prs_0831.t2.5.best', header = TRUE)
prs <- prs[,c(2,4)]
# make IDs the same in the symptom data and merge
bpm_qc <- bpm_env %>% rename('IID'='src_subject_id')
symp_gen <- merge(bpm_qc, prs, by = 'IID')
length(unique(symp_gen$IID)) # N=4135 with symptom and genetic data (EUR)

# normalise PRS to [0,1] ?

################### Impute missing data using MICE ###################

# three visualisations to inspect missing data
md.pattern(symp_gen, rotate.names = TRUE)
aggr_plot <- aggr(symp_gen, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(symp_gen), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
marginplot(symp_gen[c(1,2)])?

# impute
imputed_df <- mice(symp_gen,m=5,maxit=50,meth='pmm',seed=500)
summary(imputed_df)

completed_df <- complete(imputed_df,1)

################### Mixed graphical model for mixed data types ###################

# empty list
prs_network_list <- list()

# for each wave estimate the network and add results to list 
for (wave in c(2,4,6,8)) {
  timepoint <- completed_df$eventname == wave
  data_subset <- completed_df[timepoint, 3:ncol(completed_df)]
  prs_network <- estimateNetwork(data_subset, default="mgm")
                                 #type=c(rep('c',10), rep('g',2)),
                                 #level=c(rep(2,9), 10, rep(1,2)))
  prs_network_list[[as.character(wave)]] <- prs_network
}

## ! check warnings here

# Loop through the results and generate plot for estimated networks
for (wave in names(prs_network_list)) {
  results <- prs_network_list[[wave]]
  qgraph(results$graph, layout='spring', theme='colorblind', labels=c(unlist(labels),unlist(env_labels),"PRS"))
  title(paste("Wave =", wave), adj=0.8, line = -1.0) 
}

## try qgraph::flow with PRS? needs a qgraph object?

################### LMI testing for mplus ################### 

# longitudinal measurement invariance format (LMI) - wide 
bpm_wide <- pivot_wider(
  data = bpm_int,
  id_cols = src_subject_id,  # Identifier column
  names_from = eventname,    # Column to spread
  values_from = c(bpm_9_y, bpm_11_y, bpm_12_y, bpm_13_y, bpm_18_y, bpm_19_y)  # Values to fill the columns
)

# NAs must be -9999 for mplus
bpm_wide[is.na(bpm_wide)] <- -9999
# save
write.table(bpm_wide, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/
            Edinburgh/gmm/gmm_abcd/mplus_data/bpm_symptoms_wide.dat", col.names = FALSE)
# omit Id column
bpm_wide_3point <- bpm_wide[,c(2:49)]
# save
write.table(bpm_wide_3point, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/
            Edinburgh/networks/mplus_ri_clpm/bpm_symptoms_wide_3point.txt", col.names = FALSE)

# for RI-CLPM network (binarise symptoms)
bpm_binary <- bpm_wide %>%
  mutate(across(2:49, ~case_when(
    . == 0 ~ 0,
    . == 1 ~ 1,
    . == 2 ~ 1,
    TRUE ~ .
  )))

# checks - should be 11726 (including all missing)
length(unique(bpm_binary$src_subject_id))
# omit ID column
bpm_binary_symptoms_only <- bpm_binary[,c(2:49)]
# save
write.table(bpm_binary_symptoms_only, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/bpm_symptoms_binary.dat", col.names = FALSE, row.names = FALSE)


###### full symptoms ######

'9 I feel worthless or inferior
11 I am too fearful or anxious	
12 I feel too guilty	
13 I am self-conscious or easily embarrassed	
18 I am unhappy, sad, or depressed	
19 I worry a lot'

