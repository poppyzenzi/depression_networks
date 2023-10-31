library(dplyr)
library(qgraph)
library(tidyr)
library(bootnet)

# split into QC and network scripts 

'9 I feel worthless or inferior
11 I am too fearful or anxious	
12 I feel too guilty	
13 I am self-conscious or easily embarrassed	
18 I am unhappy, sad, or depressed	
19 I worry a lot'

labels <- c("worthless","anxious","guilty","self-conscious","unhappy","worry")

bpm <- read.csv('/Volumes/igmm/GenScotDepression/data/abcd/release5.0/core/mental-health/mh_y_bpm.csv')
symptoms = c('bpm_9_y','bpm_11_y','bpm_12_y','bpm_13_y','bpm_18_y','bpm_19_y')
int_cols = c('src_subject_id','eventname','bpm_y_scr_internal_r', symptoms)

# select cols for internalising, only for values 0-2 in symptoms (valid)
bpm_int <- bpm %>% 
  subset(select = int_cols) %>%
  filter(if_all(all_of(symptoms), ~ . >=0 &if_all(all_of(symptoms), ~ . <=2))) %>%
  mutate(`eventname` = dplyr::recode(`eventname`,
                              `baseline_year_1_arm_1`="0",
                              `6_month_follow_up_arm_1`="1",
                              `1_year_follow_up_y_arm_1`="2",
                              `18_month_follow_up_arm_1`="3",
                              `2_year_follow_up_y_arm_1`="4",
                              `30_month_follow_up_arm_1`="5",
                              `3_year_follow_up_y_arm_1`="6",
                              `42_month_follow_up_arm_1`="7",
                              `4_year_follow_up_y_arm_1`="8"))

# make sum score
#bpm_int$sumscore <- rowSums(bpm_int[symptoms])

#####################################
########## symptom tables ###########

# frequency table
freq_dat <- bpm_int
colnames(freq_dat) <- c("id","time", "sum", unlist(labels))

frequency_table <- freq_dat %>%
  gather(symptom, value, `worthless`:`worry`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time/sweep
  summarize(count = sum(value)) %>%
  spread(time, count, fill = 0)

print(frequency_table)
setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")
write.csv(frequency_table, file="abcd_symptom_frequencies.csv")

# proportion table
prop_table <- freq_dat %>%
  gather(symptom, value, `worthless`:`worry`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time
  summarize(proportion = sum(value == 1) / n()) %>% # calculate proportion of "yes" responses
  spread(time, proportion, fill = 0)

print(prop_table)
write.csv(prop_table, file="abcd_symptom_proportions.csv")

###################################
## 1: Estimate network
###################################

theme = 'Reddit'

# estimate network and loop over sweeps, use ggmModSelect or EBICglasso?
network_list <- list()

for (wave in c(2,4,6,8)) {
  timepoint <- bpm_int$eventname == wave
  data_subset <- bpm_int[timepoint, 4:ncol(bpm_int)]
  network <- estimateNetwork(data_subset,default="ggmModSelect")
  network_list[[as.character(wave)]] <- network
}

# Loop through the results and generate plot for estimated networks
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  qgraph(results$graph, layout='spring', theme=theme, labels = labels)
  #centralityPlot(results, include = c("Betweenness","Closeness","Strength"), labels=labels)
  title(paste("Wave =", wave), adj=1.0, line= -1.0) 
}

######################################
## 2: Run bootnet resampled network
#######################################

boot1_result_list <- list()

## for each wave run bootnet and store results in list 
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  boot1<-bootnet(results, nBoots=1000, default="ggmModSelect", nCores=1)
  boot1_result_list[[as.character(wave)]] <- boot1
}

## plot bootstrapped network results
for (wave in names(boot1_result_list)) {
  results <- boot1_result_list[[wave]]
  plot(results$sample, label = labels, theme=theme)
  title(paste("Wave =", wave), adj=1.0, line= -1.0) 
}

# plot bootstrapped edge CIs
# Create a list of plots and arrange in 1xn grid
plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, labels = FALSE, order = "sample")})
grid.arrange(grobs = plot_list, ncol = 4)

# bootstrapped difference tests between non-zero edge-weights
# grey = non sig diff, black = sig diff, coloured = colour of edge
## plot significant differences (alpha=0.05) of edges: 
# Create a list of edge difference plots and arrange in 1xn grid
edge_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample", labels=FALSE)})
grid.arrange(grobs = edge_plot_list, ncol = 4)

## Plot significant differences (alpha=0.05) of nodestrength
# grey = non sig, black = sig, white = value of node strength 
# Create a list of strength difference plots and arrange in 1xn grid
strength_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "strength", plot = "difference", order = "sample")})
grid.arrange(grobs = strength_plot_list, ncol = 4)

# test for difference in strength between two nodes
s1 = "bpm_19_y"
s2 = "bpm_12_y"
differenceTest(boot1, s1, s2, "strength")

##################################################
################# genetic data ###################
##################################################

setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_bpm_OUT')
prs <- read.table('t2/abcd_mdd_prs_0320.best', header = TRUE)
prs <- prs[,c(2,4)]
# make IDs the same in the symptom data and merge
bpm_qc <- bpm_int %>% rename('IID'='src_subject_id')
symp_gen <- merge(bpm_qc, prs, by = 'IID')
length(unique(symp_gen$IID)) # N=4135 with symptom and genetic data

# normalise PRS to [0,1] ?
#symp_gen$PRS <- (scale(symp_gen$PRS) - min(scale(symp_gen$PRS))) / (max(scale(symp_gen$PRS)) - min(scale(symp_gen$PRS)))

## mixed graphical model for mixed data types
prs_network_list <- list()

for (wave in c(2,4,6,8)) {
  timepoint <- symp_gen$time == wave
  data_subset <- symp_gen[timepoint, 4:ncol(symp_gen)]
  prs_network <- estimateNetwork(symp_gen[,4:ncol(symp_gen)], default="EBICglasso")
  prs_network_list[[as.character(wave)]] <- prs_network
}

# Loop through the results and generate plot for estimated networks
for (wave in names(prs_network_list)) {
  results <- prs_network_list[[wave]]
  qgraph(results$graph, layout='spring', theme='colorblind', labels=c(unlist(labels),"PRS"))
  title(paste("Wave =", wave), adj=0.8)  # line = -0.1 to lower
}


####################################
############# for mplus  ###########
####################################

# longtiudinal measurement invariance format (LMI)

bpm_wide <- pivot_wider(
  data = bpm_int,
  id_cols = src_subject_id,  # Identifier column
  names_from = eventname,    # Column to spread
  values_from = c(bpm_9_y, bpm_11_y, bpm_12_y, bpm_13_y, bpm_18_y, bpm_19_y)  # Values to fill the columns
)

bpm_wide[is.na(bpm_wide)] <- -9999

write.table(bpm_wide, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/
            Edinburgh/gmm/gmm_abcd/mplus_data/bpm_symptoms_wide.dat", col.names = FALSE)

bpm_wide_3point <- bpm_wide[,c(2:49)]

write.table(bpm_wide_3point, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/bpm_symptoms_wide_3point.txt", col.names = FALSE)


# for RI-CLPM network (make symptoms binary)
bpm_binary <- bpm_wide %>%
  mutate(across(2:49, ~case_when(
    . == 0 ~ 0,
    . == 1 ~ 1,
    . == 2 ~ 1,
    TRUE ~ .
  )))

# checks - should be 11726 (including all missing)
length(unique(bpm_binary$src_subject_id))

# drop all irrelevant cols
bpm_binary_symptoms_only <- bpm_binary[,c(2:49)]

write.table(bpm_binary_symptoms_only, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/bpm_symptoms_binary.dat", col.names = FALSE, row.names = FALSE)






