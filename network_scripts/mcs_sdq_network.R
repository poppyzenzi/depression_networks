library(dplyr)
#library(psychTools)
library(qgraph)
library(tidyr)
library(bootnet)
library(gridExtra)

setwd('/Volumes/igmm/GenScotDepression/users/poppy/mcs')
sdq_wide <- read.table('symptom_data/mcs_sdq_sym_wide.txt')
sdq_long <- read.table('symptom_data/mcs_sdq_sym_long.txt')

labels <- c("malaise", 
              "worries", 
                 "unhappy", 
                 "anxiety", 
                 "fears", 
                 "solitary",
                 "friends*",
                 "liked*",
                 "bullied",
                 "adult-oriented")
# * reverse scored
legend <- c("* reverse scored")

#### symptom frequency table #####
freq_dat <- sdq_long
colnames(freq_dat) <- c("id","sweep", unlist(labels))

frequency_table <- freq_dat %>%
  gather(symptom, value, `malaise`:`adult-oriented`) %>% # gather symptoms into key value pairs 
  group_by(sweep, symptom) %>% # group by time/sweep
  summarize(count = sum(value)) %>%
  spread(sweep, count, fill = 0)

print(frequency_table)
write.csv(frequency_table, file="symptom_data/mcs_symptom_frequencies.csv")

############################## 
######## run bootnet ########

# empty list
result_list <- list()

# run bootnet loop over sweeps
for (sweep_value in c(5,6,7)) {
  timepoint <- sdq_long$Sweep == sweep_value
  data_subset <- sdq_long[timepoint, 3:12]
  results <- bootnet(data_subset, nBoots=100, default="ggmModSelect")
  result_list[[as.character(sweep_value)]] <- results
}

# Loop through the results and generate plots
for (sweep_value in names(result_list)) {
  results <- result_list[[sweep_value]]
  plot(results$sample, label = labels)
  #title(paste("Sweep =", sweep_value), adj=0.8)
}

# ===============================

## steps for assessing accuracy of psychological network structures (Fried tutorial paper)
# a) estimate CIs on edge weights. Use bootstrapping to repeatedly estimate model under sampled data and estimate statistic of interest.
# b) assess stability of centrality indices in subsets of cases
# c) test for sig diff btw edge weights and centrality indices

# symptoms are ordinal: compute polychoric correlation matrix as input using qgraph::cor_auto

network_list <- list()

# estimate network and loop over sweeps, use ggmModSelect or EBICglasso?
for (sweep_value in c(5,6,7)) {
  timepoint <- sdq_long$Sweep == sweep_value
  data_subset <- sdq_long[timepoint, 3:12]
  network <- estimateNetwork(data_subset,default="ggmModSelect")
  network_list[[as.character(sweep_value)]] <- network
}

theme = 'Reddit'
# Loop through the results and generate plot for estimated networks
for (sweep_value in names(network_list)) {
  results <- network_list[[sweep_value]]
  #qgraph(results$graph, layout='spring', theme=theme, labels = labels)
  centralityPlot(results, include = c("Betweenness","Closeness","Strength"), labels=labels)
  title(paste("Sweep =", sweep_value), adj=0.8)  # line = -0.1 to lower
}

boot1_result_list <- list()

## non-parametric bootstrap (important to use for ordinal data)
# 1000 boots quite slow
for (sweep_value in names(network_list)) {
  results <- network_list[[sweep_value]]
  boot1<-bootnet(results, nBoots=1000, default="ggmModSelect", nCores=1)
  # plot bootstrapped edgeCIs: indicates care to be taken interpreting edges in network 
  boot1_result_list[[as.character(sweep_value)]] <- boot1
  
}

## plot bootstrapped edge CIs
# Create a list of plots and arrange in 1x3 grid
plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, labels = FALSE, order = "sample")})
# Arrange the plots in a 1x3 grid
grid.arrange(grobs = plot_list, ncol = 3)


## plot significant differences (alpha=0.05) of edges: 
# Create a list of edge difference plots and arrange in 1x3 grid
edge_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample", labels=FALSE)})
grid.arrange(grobs = edge_plot_list, ncol = 3)


## Plot significant differences (alpha=0.05) of nodestrength:
# Create a list of strength difference plots and arrange in 1x3 grid
strength_plot_list <- lapply(boot1_result_list, function(boot_result) {
  plot(boot_result, "strength", plot = "difference", order = "sample")})
grid.arrange(grobs = strength_plot_list, ncol = 3)


# test for difference in strength between node x and y
differenceTest(boot1,"PSDHS00","PSDNC00","strength")

