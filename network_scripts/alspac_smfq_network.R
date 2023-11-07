### Rscript for making symptom and risk factor networks in ALPSAC

library(haven)
library(tidyr)
library(dplyr)
library(purrr)
library(data.table)
library(bootnet)
library(psychonetrics)
library(gridExtra)
library(qgraph)

## load symptom data (from Rscript: alspac_smfq_QC.R)

setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
smfq_dat <- read.table('smfq_sypmtoms_wide.txt', check.names = FALSE)

# check right order of symptom labels
labels <- c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
            "tearful", "distracted", "self-loathing", "guilty", "isolated", 
            "unloved", "inadequate", "incompetent")


############# quality control ###############

## filter for only valid symptom scores btw 0-3
smfq_symptoms <- smfq_dat %>%
  filter(if_all(.cols = all_of(names(reshaped)[3:ncol(smfq_symptoms)]), ~ . >= 0 & . <= 3))

## SMFQ: 1 = true, 2 = sometimes, 3 = not
## recode scores to binary 1 true 0 not true
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

# save
write.table(smfq_symptoms, file='smfq_symptoms_qcd.txt')


########## symptom tables ###########

# frequency table
freq_dat <- smfq_symptoms
colnames(freq_dat) <- c("id","time", unlist(labels))

frequency_table <- freq_dat %>%
  gather(symptom, value, `unhappy`:`incompetent`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time/sweep
  summarize(count = sum(value)) %>%
  spread(time, count, fill = 0)

print(frequency_table)
write.csv(frequency_table, file="alspac_symptom_frequencies.csv")

# proportion/endorsement table
prop_table <- freq_dat %>%
  gather(symptom, value, `unhappy`:`incompetent`) %>% 
  group_by(time, symptom) %>% # group by time/sweep
  summarize(proportion = sum(value == 1) / n()) %>% # calculate proportion of "yes" responses
  spread(time, proportion, fill = 0)

print(prop_table)
write.csv(prop_table, file="alspac_symptom_proportions.csv")

###################################
## 1: Estimate symptom only network
###################################

theme = 'Reddit'

# estimate network and loop over sweeps, use ggmModSelect or EBICglasso?
network_list <- list()

for (wave in c(1,3,4)) {
  timepoint <- smfq_symptoms$time == wave
  data_subset <- smfq_symptoms[timepoint, 3:ncol(smfq_symptoms)]
  network <- estimateNetwork(data_subset,default="ggmModSelect")
  network_list[[as.character(wave)]] <- network
}

# Loop through the results and generate plot for estimated networks
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  qgraph(results$graph, layout='spring', theme=theme, labels = labels)
  title(paste("Wave =", wave), adj=0.8, line=-1)
  centralityPlot(results, include = c("Betweenness","Closeness","Strength"), labels=labels)
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
  plot(results$sample, label = labels, theme=theme)
  title(paste("Wave =", wave), adj=0.8, line=-1.0)
}

# plot bootstrapped edge CIs
# Create a list of plots and arrange in 1xn grid
plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, labels = FALSE, order = "sample")})
grid.arrange(grobs = plot_list, ncol = length(boot1_result_list))


# bootstrapped difference tests between non-zero edge-weights
# grey = non sig diff, black = sig diff, coloured = colour of edge
## plot significant differences (alpha=0.05) of edges: 
# Create a list of edge difference plots and arrange in 1xn grid
edge_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample", labels=FALSE)})
grid.arrange(grobs = edge_plot_list, ncol = length(boot1_result_list))

## Plot significant differences (alpha=0.05) of nodestrength
# grey = non sig, black = sig, white = value of node strength 
# Create a list of strength difference plots and arrange in 1xn grid
strength_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "strength", plot = "difference", order = "sample")})
grid.arrange(grobs = strength_plot_list, ncol = length(boot1_result_list))

# test for difference in strength between two nodes
s1 = "01"
s2 = "10"
differenceTest(boot1, s1, s2, "strength")


################# environmental data ###################

# read in environmental vars
predictors <- read.table('alspac_envi_vars')

# merge with symptom data by id 
smfq_env <- merge(smfq_symptoms, predictors, by='id')

smfq_env <- smfq_env %>% rename(`sex` = kz021,
                                `maternal depression` = r2021,
                                `bullying` = f8fp470,
                                `child trauma` = AT5_n,
                                `sleep` = FJCI250,
                                `income` = h470,
)

env_labels = list(names(smfq_env[,16:ncol(smfq_env)]))

################# genetic data ###################

setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/prs/prs_alspac_OUT')
prs <- read.table('alspac_mdd_prs_0320.best', header = TRUE)
prs <- prs[,c(2,4)]

prs <- prs %>% rename(`MDD` = `PRS`)
                                
# make IDs the same in the symptom data and merge
smfq_qc <- smfq_env %>% mutate(id = gsub("_", "", id)) %>% rename('IID'='id')
symp_gen <- merge(smfq_qc, prs, by = 'IID')
length(unique(symp_gen$IID)) # N=6096 with symptom and genetic data

# normalise PRS to mean of 0 and SD of 1 
symp_gen$MDD <- scale(symp_gen$MDD, center = TRUE, scale = TRUE)

######################################################################
###### run mixed graphical model network for mixed data types ########
###################################################################### 

all_network_list <- list()

for (wave in c(1,3,4)) {
  timepoint <- symp_gen$time == wave
  data_subset <- symp_gen[timepoint, 3:ncol(symp_gen)]
  all_network <- estimateNetwork(data_subset, default="EBICglasso")
  all_network_list[[as.character(wave)]] <- all_network
}

# network colours: can be "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
theme <- 'colorblind'
# node colours: can be "rainbow","colorblind", "pastel","gray","R","ggplot2"
palette <- 'pastel'

groups=list(1:13,14:19,20)
group_names <- c("SMFQ symptoms", "Environmental factors", "Polygenic Risk Scores")
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

