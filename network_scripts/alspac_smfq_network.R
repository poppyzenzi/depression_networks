library(haven)
library(tidyr)
library(dplyr)
library(purrr)
library(data.table)
library(bootnet)
library(psychonetrics)
library(gridExtra)


#########################

# load symptom data
# from Rscript: alspac_smfq_QC.R

setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
smfq_dat <- read.table('smfq_sypmtoms_wide.txt', check.names = FALSE)

# SMFQ: 1 = true, 2 = sometimes, 3 = not

# check right order
labels <- c("unhappy", 
            "anhedonia", 
            "apathetic", 
            "restless", 
            "worthless", 
            "tearful",
            "distracted",
            "self-loathing",
            "guilty",
            "isolated",
            "unloved",
            "inadequate",
            "incompetent")


############# quality control ###############

## filter for only valid symptom scores btw 0-3
smfq_symptoms <- smfq_dat %>%
  filter(if_all(.cols = all_of(names(reshaped)[3:ncol(smfq_symptoms)]), ~ . >= 0 & . <= 3))

## decide here if want to make binary 
smfq_binary <- smfq_symptoms %>%
  mutate(across(3:15, ~case_when(
    . == 1 ~ 1,
    . == 2 ~ 1,
    . == 3 ~ 0,
    TRUE ~ .
  )))

smfq_symptoms <- smfq_binary

# make time col numeric
mapping <- c("t1" = 1, "t2" = 2, "t3" = 3, "t4" = 4)
smfq_symptoms$time <- mapping[smfq_symptoms$time]

#####################################

#### symptom frequency table #####
freq_dat <- smfq_symptoms
colnames(freq_dat) <- c("id","time", unlist(labels))

frequency_table <- freq_dat %>%
  gather(symptom, value, `unhappy`:`incompetent`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time/sweep
  summarize(count = sum(value)) %>%
  spread(time, count, fill = 0)

print(frequency_table)
write.csv(frequency_table, file="alspac_symptom_frequencies.csv")

prop_table <- freq_dat %>%
  gather(symptom, value, `unhappy`:`incompetent`) %>% # explicitly specify the order of symptoms
  group_by(time, symptom) %>% # group by time/sweep
  summarize(proportion = sum(value == 1) / n()) %>% # calculate proportion of "yes" responses
  spread(time, proportion, fill = 0)

print(prop_table)
write.csv(prop_table, file="alspac_symptom_proportions.csv")


####### run bootnet sample ########

# empty list
result_list <- list()

# run bootnet loop over sweeps
for (wave in c(1,2,3,4)) {
  timepoint <- smfq_symptoms$time == wave
  data_subset <- smfq_symptoms[timepoint, 3:ncol(smfq_symptoms)]
  results <- bootnet(data_subset, nBoots=20, default="ggmModSelect")
  result_list[[as.character(wave)]] <- results
}

# Loop through the results and generate plots
for (wave in names(result_list)) {
  results <- result_list[[wave]]
  plot(results$sample, label = labels)
  title(paste("Wave =", wave), adj=0.8)  # line = -0.1 to lower
}

######################################################

## steps for assessing accuracy of psychological network structures (Fried tutorial paper)
# a) estimate CIs on edge weights. Use bootstrapping to repeatedly estimate model under sampled data and estimate statistic of interest.
# b) assess stability of centrality indices in subsets of cases
# c) test for sig diff btw edge weights and centrality indices

# symptoms are ordinal: compute polychoric correlation matrix as input using qgraph::cor_auto

# estimate network and loop over sweeps, use ggmModSelect or EBICglasso?
network_list <- list()

for (wave in c(1,2,3,4)) {
  timepoint <- smfq_symptoms$time == wave
  data_subset <- smfq_symptoms[timepoint, 3:ncol(smfq_symptoms)]
  network <- estimateNetwork(data_subset,default="ggmModSelect")
  network_list[[as.character(wave)]] <- network
}

theme = 'Reddit'

# Loop through the results and generate plot for estimated networks
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  qgraph(results$graph, layout='spring', theme=theme, labels = labels)
  #centralityPlot(results, include = c("Betweenness","Closeness","Strength"), labels=labels)
  title(paste("Wave =", wave), adj=0.8)  # line = -0.1 to lower
}

## non-parametric bootstrap (important to use for ordinal data)
# 1000 boots quite slow
boot1_result_list <- list()

for (wave in names(network_list)) {
  results <- network_list[[wave]]
  boot1<-bootnet(results, nBoots=100, default="ggmModSelect", nCores=1)
  # plot bootstrapped edgeCIs: indicates care to be taken interpreting edges in network 
  boot1_result_list[[as.character(wave)]] <- boot1
  
}

## plot bootstrapped edge CIs
# Create a list of plots and arrange in 1x3 grid
plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, labels = FALSE, order = "sample")})
# Arrange the plots in a 1x4 grid
grid.arrange(grobs = plot_list, ncol = 4)


## plot significant differences (alpha=0.05) of edges: 
# Create a list of edge difference plots and arrange in 1x3 grid
edge_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample", labels=FALSE)})
grid.arrange(grobs = edge_plot_list, ncol = 4)


## Plot significant differences (alpha=0.05) of nodestrength:
# Create a list of strength difference plots and arrange in 1x3 grid
strength_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "strength", plot = "difference", order = "sample")})
grid.arrange(grobs = strength_plot_list, ncol = 4)


# test for difference in strength between two nodes
s1 = ""
s2 = ""
differenceTest(boot1, s1, s2, "strength")

################################################################
####################### genetic data ###########################
################################################################

labels <- c("unhappy", 
            "anhedonia", 
            "apathetic", 
            "restless", 
            "worthless", 
            "tearful",
            "distracted",
            "self-loathing",
            "guilty",
            "isolated",
            "unloved",
            "inadequate",
            "incompetent",
            "PRS")

setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_alspac_OUT')
prs <- read.table('alspac_mdd_prs_0320.best', header = TRUE)
prs <- prs[,c(2,4)]

# make IDs the same in the symptom data
smfq_qc <- smfq_symptoms %>% mutate(id = gsub("_", "", id)) %>% rename('IID'='id')
# merge on IID
symp_gen <- merge(smfq_qc, prs, by = 'IID')
length(unique(symp_gen$IID)) # N=6096 with symptom and genetic data
                
## fix scoring 
# currently: 1 = true, 2 = sometimes, 3 = not
# can also make binary here
symp_gen <- symp_gen %>%
   mutate(across(3:15, ~case_when(
    . == 1 ~ 1,
    . == 2 ~ 1,
    . == 3 ~ 0,
    TRUE ~ .
  )))

# normalise PRS to [0,1]
#symp_gen$PRS <- (scale(symp_gen$PRS) - min(scale(symp_gen$PRS))) / (max(scale(symp_gen$PRS)) - min(scale(symp_gen$PRS)))

## mixed graphical model for mixed data types

prs_network <- estimateNetwork(symp_gen[,3:ncol(symp_gen)], default="EBICglasso")
qgraph(prs_network$graph, layout='spring', theme='colorblind', labels = labels)





###### run bootnet for sample ######

# empty list
result_list <- list()

# run bootnet loop over sweeps
for (wave in c(1,2,3,4)) {
  timepoint <- symp_gen$time == wave
  data_subset <- symp_gen[timepoint, 3:ncol(symp_gen)]
  results <- bootnet(data_subset, nBoots=20, default="ggmModSelect", covar)
  result_list[[as.character(wave)]] <- results
}

# Loop through the results and generate plots
for (wave in names(result_list)) {
  results <- result_list[[wave]]
  plot(results$sample, label = labels)
  title(paste("Wave =", wave), adj=0.8)  # line = -0.1 to lower
}


## estimate network 

network_list <- list()

for (wave in c(1,2,3,4)) {
  timepoint <- symp_gen$time == wave
  data_subset <- symp_gen[timepoint, 3:ncol(symp_gen)]
  network <- estimateNetwork(data_subset,default="ggmModSelect")
  network_list[[as.character(wave)]] <- network
}

theme = 'Reddit'

# Loop through the results and generate plot for estimated networks
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  qgraph(results$graph, layout='spring', pastel=TRUE, labels = labels)
  #centralityPlot(results, include = c("Betweenness","Closeness","Strength"), labels=labels)
  title(paste("Wave =", wave), adj=0.8)  # line = -0.1 to lower
}

