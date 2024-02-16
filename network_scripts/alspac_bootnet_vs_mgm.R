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
library(viridis)
library(mgm)

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

## or keep trichotomous
smfq_tri <- smfq_symptoms %>%
  mutate(across(3:15, ~case_when(
    . == 1 ~ 2,
    . == 2 ~ 1,
    . == 3 ~ 0,
    TRUE ~ .
  )))

smfq_symptoms <- smfq_tri

# make time col numeric
mapping <- c("t1" = 1, "t2" = 2, "t3" = 3, "t4" = 4)
smfq_symptoms$time <- mapping[smfq_symptoms$time]

# save
write.table(smfq_symptoms, file='smfq_symptoms_qcd.txt')

################### ################### ################
################# environmental data ###################

# read in environmental vars
setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
predictors <- read.table('alspac_envi_vars.txt')

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

# save df with all vars
write.table(symp_gen, '/Volumes/igmm/GenScotDepression/users/poppy/alspac/network_all_vars.txt')


##### add class data

# Set the working directory
setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_alspac')
# Read data from a text file
alspac_4k <- read.table('mplus_data/4k_gmm_smfq_gen_only.txt', header = FALSE)
# Rename columns
colnames(alspac_4k) <- c('y0', 'y1', 'y2', 'y3', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6', 'v7', 'v8', 'v9', 'v10', 'class', 'SubjectNumeric')
alspac_4k <- alspac_4k[, c('SubjectNumeric', 'class')]
alspac_4k[alspac_4k == '*'] <- NA
# merge with IID
alspac_ids <- read.csv("genetic_subjects_only/alspac_dep_long_gen_only.csv")
alspac_ids <- alspac_ids[, c('IID', 'SubjectNumeric')]
# Merge the 'df' and 'alspac_4k' data frames on the 'SubjectNumeric' column
class_data <- merge(alspac_ids, alspac_4k, by ="SubjectNumeric")
class_data <- class_data[, c('IID', 'class')]
# drop duplicate rows as 4 per IID
class_data <- unique(class_data)

# finally, merge class data with network dataframe
symp_gen_class <- merge(symp_gen, class_data, by='IID')

symp_gen_class <- symp_gen_class %>%
  select(IID, time, class, everything())

symp_gen_class <- symp_gen_class %>% 
  mutate(class = case_when(
    class == 1 ~ "Decreasing",
    class == 2 ~ "Stable low",
    class == 3 ~ "Persistent",
    class == 4 ~ "Increasing",
    TRUE ~ as.character(class)  # Keep other values as is
  ))

colnames(symp_gen_class) <- c('IID', 'time', 'class', 1:20)

complete_rows <- complete.cases(symp_gen_class[, 4:16])
symp_gen_class_complete <- symp_gen_class[complete_rows, ]

symp_gen_class_complete_filtered <- symp_gen_class_complete %>%
  group_by(IID) %>%
  filter(all(c('1', '2', '3', '4') %in% time)) %>%
  ungroup()

symp_gen_complete <- symp_gen %>%
  group_by(IID) %>%
  filter(all(c('1', '2', '3', '4') %in% time)) %>%
  ungroup()

##########################################################################
## bootnet::estimateNetwork with mixed data types using graphical lasso ##
########################################################################## 

# network colours: can be "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
theme <- 'colorblind'
# node colours: can be "rainbow","colorblind", "pastel","gray","R","ggplot2"
palette <- 'pastel'
# group quality control 
groups=list(1:13,14:19,20)
group_names <- c("SMFQ symptoms", "Environmental factors", "Genetic factors")
names(groups) <- group_names

all_network_list <- list()

# first estimate the network using bootnet::estimateNetwork
set.seed(2714658)
for (wave in c(1,2,3,4)) {
  timepoint <- symp_gen_complete$time == wave
  data_subset <- symp_gen_complete[timepoint, 3:ncol(symp_gen_complete)]
  all_network <- estimateNetwork(data_subset, default="mgm",
                                 type = c(rep('c', 19), 'g'), # c = cat, g = gaussian [check with hist()]
                                 level = c(rep(3, 13), rep(2, 4), rep(5, 2), 1)) # levels, 1 for continuous
  all_network_list[[as.character(wave)]] <- all_network
}

# stratify networks by latent class
class_values <- c("Decreasing", "Stable low", "Persistent", "Increasing")
for (wave in c(1, 2, 3, 4)) {
  for (class_value in class_values) {
    timepoint <- symp_gen_class_complete_filtered$time == wave
    class_condition <- symp_gen_class_complete_filtered$class == class_value
    data_subset <- symp_gen_class_complete_filtered[timepoint & class_condition, 4:ncol(symp_gen_class_complete_filtered)]
    all_network <- estimateNetwork(data_subset, default="mgm",
                                   type = c(rep('c', 19), 'g'),
                                   level = c(rep(3, 13), rep(2, 4), rep(5, 2), 1))
    all_network_list[[paste0("Wave", wave, "_Class", class_value)]] <- all_network
  }
}

# Loop through the results and generate plot for estimated networks
for (wave in names(all_network_list)) {
  results <- all_network_list[[wave]]
  qgraph(results$graph, 
         layout='spring', repulsion = 1, 
         theme=theme, groups=groups, palette = palette, 
         nodeNames=c(unlist(labels),c(unlist(env_labels), 'PRS')),
         legend=TRUE, legend.mode='style1', legend.cex=0.4,
        vsize = 5.0, esize = 15)
  title(paste('Wave = ', wave), adj=0, line=-1, cex.main=1.2)
  #qgraph::centralityPlot(results, include = c("Strength","Closeness","Betweenness"), 
  #labels=c(unlist(labels),c(unlist(env_labels), 'PRS')),
  #orderBy = 'Strength', theme_bw = TRUE)
}


#####

# Set the number of rows and columns for the plot grid
n_rows <- ceiling(length(all_network_list) / 4)  
n_cols <- 4  # Number of columns

# Create a multi-plot layout
par(mfrow = c(n_rows, n_cols))

# Loop through the results and generate individual network plots
for (wave in names(all_network_list)) {
  results <- all_network_list[[wave]]
  qgraph(results$graph, 
         layout='spring', repulsion = 1, 
         theme=theme, groups=groups, palette = palette, 
         nodeNames=c(unlist(labels),c(unlist(env_labels), 'PRS')),
         legend=FALSE, legend.mode='style1', legend.cex=0.4,
         vsize = 5.0, esize = 15)
  title(paste(wave), line=0, cex.main=0.8)  # Set adj to 0, line to -2, and reduce cex.main
}

# Reset the plotting layout to its default value
par(mfrow = c(1, 1))




#####


##### Assessing stability: ##################
##### 2: Run bootnet resampled network ######

all_boot_result_list <- list()

## for each wave run bootnet and store results in list 
for (wave in names(all_network_list)) {
  results <- all_network_list[[wave]]
  boot1<-bootnet(results, nBoots=20, default="EBICglasso", nCores=1)
  all_boot_result_list[[as.character(wave)]] <- boot1
}

## plot bootstrapped network results
for (wave in names(all_boot_result_list)) {
  results <- all_boot_result_list[[wave]]
  plot(results$sample, label = c(unlist(labels),c(unlist(env_labels),'PRS')), theme=theme)
  title(paste("Wave =", wave), adj=0.8, line=-1.0)
}

# plot bootstrapped edge CIs
# Create a list of plots and arrange in 1xn grid
plot_list <- lapply(all_boot_result_list, function(boot_result) {
  plot(boot_result, labels = FALSE, order = "sample")})
grid.arrange(grobs = plot_list, ncol = length(all_boot_result_list))

# bootstrapped difference tests between non-zero edge-weights
# grey = non sig diff, black = sig diff, coloured = colour of edge
## plot significant differences (alpha=0.05) of edges: 
# Create a list of edge difference plots and arrange in 1xn grid
edge_plot_list <- lapply(all_boot_result_list, function(boot_result) {
  plot(boot_result, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample", labels=FALSE)})
grid.arrange(grobs = edge_plot_list, nrow=2, ncol=2)

## Plot significant differences (alpha=0.05) of nodestrength
# grey = non sig, black = sig, white = value of node strength 
# Create a list of strength difference plots and arrange in 1xn grid
strength_plot_list <- lapply(all_boot_result_list, function(boot_result) {
  plot(boot_result, "strength", plot = "difference", order = "sample",
       labels=TRUE, legend=TRUE
  )})
grid.arrange(grobs = strength_plot_list, nrow=2, ncol=2)

########################################################
#### mgm for mixed data types using graphical lasso ####
########################################################

# group QC for plotting
groups=list(1:13,14:19,20)
group_names <- c("SMFQ symptoms", "Environmental factors", "Genetic factors")
names(groups) <- group_names
categorical_palette <- viridis_pal(option = "D")(8)
group_colours <- c(categorical_palette[1],
                   categorical_palette[4],
                   categorical_palette[8])

colnames(symp_gen)[3:15] <- labels
data <- na.omit(symp_gen_complete[,3:ncol(symp_gen_complete)])

# fit model at 4 waves
fit_list <- list()
for (wave in c(1,2,3,4)) {
  timepoint <- symp_gen_complete$time == wave
  data_subset <- na.omit(symp_gen_complete[timepoint, 3:ncol(symp_gen_complete)])
  fit_alspac <- mgm(data = as.matrix(data_subset), 
                 type = c(rep('c', 19), 'g'),
                 level = c(rep(2, 17), rep(5, 2), 1),
                 k = 2, 
                 lambdaSel = 'EBIC', 
                 lambdaGam = 0.25)
  fit_list[[as.character(wave)]] <- fit_alspac
}

# grey edges = for cat vars where no sign is defined
# edge width proportional to abs value of edge-parameter
for (wave in names(fit_list)) {
    results <- fit_list[[wave]]$pairwise$wadj
    qgraph(results, 
           layout = 'spring', repulsion = 1,
           edge.color = fit_alspac$pairwise$edgecolor, 
           nodeNames = names(data),
           groups = groups,
           legend.mode="style2", legend.cex=0.42, 
           vsize = 3.5, esize = 15,
           theme=theme,
           palette = palette, legend.mode='style1')
    title(paste("Wave = ", wave), adj=0.1, line=-1.5)
}
