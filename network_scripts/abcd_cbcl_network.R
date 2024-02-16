### network for CBCL DSM-5 in ABCD

library(dplyr)
library(qgraph)
library(tidyr)
library(bootnet)

setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")

cbcl_binary <- read.table('cbcl_symptoms_binary.txt')

labels <- c("anhedonia","cries","harms self","poor eating","worthless",
            "guilty", "tired", "sleeps less", "sleeps more", "suicidal", 
            "sleep problems", "lacks energy", "sad")


#####################################
########## symptom tables ###########

# frequency table
freq_dat <- cbcl_binary
colnames(freq_dat) <- c("id","time", unlist(labels))

frequency_table <- freq_dat %>%
  gather(symptom, value, `anhedonia`:`sad`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time/sweep
  summarize(count = sum(value)) %>%
  spread(time, count, fill = 0)

print(frequency_table)
write.csv(frequency_table, file="abcd_cbcl_symptom_frequencies.csv")

# proportion table
prop_table <- freq_dat %>%
  gather(symptom, value, `anhedonia`:`sad`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time
  summarize(proportion = sum(value == 1) / n()) %>% # calculate proportion of "yes" responses
  spread(time, proportion, fill = 0)

print(prop_table)
write.csv(prop_table, file="abcd_cbcl_symptom_proportions.csv")


###################################
## 1: Estimate network
###################################

theme = 'colorblind'

# estimate network and loop over sweeps, use ggmModSelect or EBICglasso?
network_list <- list()

# ages 11,12,13,14
for (wave in c(2,4,6,8)) {
  timepoint <- cbcl_binary$eventname == wave
  data_subset <- cbcl_binary[timepoint, 3:ncol(cbcl_binary)]
  network <- estimateNetwork(data_subset,default="ggmModSelect")
  network_list[[as.character(wave)]] <- network
}

# Loop through the results and generate plot for estimated networks
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  qgraph(results$graph, layout='spring', theme=theme, labels = labels)
  #centralityPlot(results, include = c("Betweenness","Closeness","Strength"), labels=labels)
  title(paste("Wave =", wave), adj=0.8, line= -1.0) 
}

######################################
## 2: Run bootnet resampled network
#######################################

boot1_result_list <- list()

## for each wave run bootnet and store results in list 
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  boot1<-bootnet(results, nBoots=100, default="ggmModSelect", nCores=1)
  boot1_result_list[[as.character(wave)]] <- boot1
}

## plot bootstrapped network results
for (wave in names(boot1_result_list)) {
  results <- boot1_result_list[[wave]]
  plot(results$sample, label = labels, theme='colorblind')
  #title(paste("Wave =", wave), adj=1.0, line= -1.0) 
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
s1 = "cbcl_q05_p"
s2 = "cbcl_q76_p"
differenceTest(boot1, s1, s2, "strength")

################ environmental data ###################

# read in environmental vars

### check time labelling of predictors ##### from 1234 to 2468
predictors <- read.table('/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_envi_vars.txt')

# recode time points
predictors_qc <- predictors %>%
  mutate(eventname = case_when(
    eventname == 1 ~ 2,
    eventname == 2 ~ 4,
    eventname == 3 ~ 6,
    eventname == 4 ~ 8,
    TRUE ~ eventname
  ))

# merge with symptom data by id 
cbcl_env <- merge(cbcl_binary, predictors_qc, by=c('src_subject_id','eventname'))

cbcl_env <- cbcl_env %>% rename(`sex` = demo_sex_v2,
                              `maternal depression` = fam_history_q6d_depression,
                              `bullying` = kbi_p_c_bully,
                              #`child trauma` = ,
                              `sleep` = sds_p_ss_total,
                              `income` = demo_comb_income_v2,
)

env_labels = list(names(cbcl_env[,16:ncol(cbcl_env)]))

##################################################
################# genetic data ###################
##################################################

setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_bpm_OUT')
prs <- read.table('t2/abcd_mdd_prs_0320.best', header = TRUE)
prs <- prs[,c(2,4)]
# make IDs the same in the symptom data and merge
cbcl_qc <- cbcl_env %>% rename('IID'='src_subject_id')
symp_gen <- merge(cbcl_qc, prs, by = 'IID')
length(unique(symp_gen$IID)) # N=4156 with symptom and genetic data

# normalise PRS to mean of 0 and SD of 1 
symp_gen$PRS <- scale(symp_gen$PRS, center = TRUE, scale = TRUE)

# save df with all vars
write.table(symp_gen, '/Volumes/igmm/GenScotDepression/users/poppy/abcd/network_cbcl_all_vars.txt')

######################################################################
###### run mixed graphical model network for mixed data types ########
###################################################################### 

all_network_list <- list()

for (wave in c(4,6,8)) {
  timepoint <- symp_gen$eventname == wave
  data_subset <- symp_gen[timepoint, 3:ncol(symp_gen)]
  all_network <- estimateNetwork(data_subset, default="EBICglasso")
  all_network_list[[as.character(wave)]] <- all_network
}

# network colours: can be "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
theme <- 'colorblind'
# node colours: can be "rainbow","colorblind", "pastel","gray","R","ggplot2"
palette <- 'pastel'

groups=list(1:13,14:19,20)
group_names <- c("SMFQ symptoms", "Environmental factors", "Genetic factors")
names(groups) <- group_names

# Loop through the results and generate plot for estimated networks
for (wave in names(all_network_list)) {
  results <- all_network_list[[wave]]
  qgraph(results$graph, layout='spring', theme=theme, groups=groups,
         palette = palette, nodeNames=c(unlist(labels),c(unlist(env_labels),'PRS')),
         legend=TRUE, legend.mode='style1', legend.cex=0.45)
  title(paste("Wave =", wave), adj=0.1, line=-1.5)
  
}

#####################################



