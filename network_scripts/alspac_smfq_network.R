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
library(mice)
library(ggmice)
library(NetworkComparisonTest)
library(viridis)

## load symptom data (from Rscript: alspac_smfq_QC.R)

setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
smfq_dat <- read.table('smfq_symptoms_wide.txt', check.names = FALSE)

# check right order of symptom labels
labels <- c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
            "tearful", "distracted", "self_loathing", "guilty", "isolated", 
            "unloved", "inadequate", "incompetent")

colnames(smfq_dat) <- c('id', 'time', unlist(labels))


############# quality control ###############

## filter for only valid symptom scores btw 0-3
smfq_symptoms <- smfq_dat %>%
  filter(if_all(.cols = all_of(names(smfq_dat)[3:ncol(smfq_dat)]), ~ . >= 0 & . <= 3))

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
mapping <- c("t1" = 1, "t2" = 2, "t3" = 3, "t4" = 4, "t5" = 5, "t6" = 6, "t7" = 7)
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

print(prop_table %>% select(1,2,4,6))

write.csv(prop_table, file="alspac_symptom_proportions.csv")
write.csv(prop_table %>% select(1,2,4,6), file="alspac_symptom_proportions_3wav.csv")

###################################
## 1: Estimate symptom only network
###################################

theme = 'Reddit'

# estimate network and loop over sweeps, use ggmModSelect or EBICglasso?
network_list <- list()

for (wave in c(1,3,5)) {
  timepoint <- smfq_symptoms$time == wave
  data_subset <- smfq_symptoms[timepoint, 3:ncol(smfq_symptoms)]
  network <- estimateNetwork(data_subset,default="ggmModSelect")
  network_list[[as.character(wave)]] <- network
}

# Loop through the results and generate plot for estimated networks
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  qgraph(results$graph, layout='spring', theme=theme, labels = labels)
  title(paste("Wave =", wave), adj=0.8, line=-0.1)
  centralityPlot(results, include = c("Betweenness","Closeness","Strength"), labels=labels)
}

######################################
## 2: Run bootnet resampled network
#######################################

boot1_result_list <- list()

## for each wave run bootnet and store results in list 
for (wave in names(network_list)) {
  results <- network_list[[wave]]
  boot1<-bootnet(results, nBoots=50, default="ggmModSelect", nCores=1)
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

################### ################### ################
################# environmental data ###################

# read in environmental vars
setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
predictors <- read.table('alspac_envi_vars.txt')

# merge with symptom data by id 
smfq_env <- merge(smfq_symptoms, predictors, by='id')

smfq_env <- smfq_env %>% rename(`sex` = kz021,
                                `maternal_depression` = r2021,
                                `bullying` = f8fp470,
                                `child_trauma` = AT5_n,
                                `sleep` = FJCI250,
                                `income` = h470,
)

env_labels = list(names(smfq_env[,16:ncol(smfq_env)]))

################# genetic data ###################

setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/projects/prs/prs_alspac_OUT/all_thresholds')
#pfac_prs <- read.table('alspac_high_prs_0817.t2.best',header = TRUE)
#pfac_prs <- pfac_prs[,c(2,4)]
#pfac_prs <- pfac_prs %>% rename(`p-factor` = `PRS`)

# reading in parcel PRS data
csvs <- c('alspac_mood_prs1221.t3.best', 'alspac_psychotic_prs1221.t3.best', 'alspac_neurodev_prs1221.t3.best')

# Create an empty list to store dataframes
prs_data_list <- list()

# iterate over each file
for (csv_file in csvs) {
  df <- read.table(csv_file, header = TRUE)
  df <- df[, c('IID', 'PRS')]
  col_prefix <- strsplit(csv_file, "_")[[1]][2]  # extract the prefix of the column name from the file name
  col_prefix2 <- toupper(col_prefix) # capitalise
  setnames(df, c('PRS'), paste0(col_prefix2))  # rename the second column to the appropriate prefix
  prs_data_list[[col_prefix]] <- df  # Store each PRS dataframe in the list
}

# Combine all PRS dataframes into a single dataframe based on IID column
combined_prs_data <- Reduce(function(x, y) merge(x, y, by = 'IID', all.x = TRUE), prs_data_list)

# make IDs the same in the symptom data and merge
smfq_qc <- smfq_env %>% mutate(id = gsub("_", "", id)) %>% rename('IID'='id')
symp_gen <- merge(smfq_qc, combined_prs_data, by = 'IID')
length(unique(symp_gen$IID)) # N=6096 with symptom and genetic data (t4) (t7 = 6240)

# normalise PRS to mean of 0 and SD of 1
cols_to_standardise <- names(symp_gen)[(ncol(symp_gen) - 2):ncol(symp_gen)]
symp_gen[, cols_to_standardise] <- scale(symp_gen[, cols_to_standardise])
# check mean and SD
summary(symp_gen[, cols_to_standardise])

# save df with all vars
write.table(symp_gen, '/Volumes/igmm/GenScotDepression/users/poppy/alspac/network_all_vars.txt')


##### add class data (optional, if want to stratify by trajectory)
'''
# Set the working directory
setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/projects/gmm/gmm_alspac')
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
'''
######################################################################
################### Impute missing data using MICE ###################
######################################################################

### Prepare the df and merge with og df
template <- expand.grid(IID = unique(symp_gen$IID), time = c(1,3,5))
prepped_df <- merge(template, symp_gen, by = c("IID", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(IID, time)
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(unhappy:income, ~replace(., is.na(.), NA)))
# fill PRS for time invariant columns
prepped_df <- prepped_df %>%
  group_by(IID) %>%
  fill(MOOD:NEURODEV, .direction = "down") %>%
  fill(MOOD:NEURODEV, .direction = "up")
# reorder cols
prepped_df <- prepped_df[, c("IID", "time", names(symp_gen)[3:ncol(symp_gen)])]
# Check completed df
print(prepped_df)

# remove individuals without genetic data
complete_rows <- complete.cases(prepped_df[, 22:24])
prepped_df <- prepped_df[complete_rows, ]

# visualisations to inspect missing data
md.pattern(prepped_df, rotate.names = TRUE)
plot_pattern(prepped_df, square = TRUE, rotate = TRUE)
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

# make predictor matrix
# we only want to impute the symptom and environmental vars
predMat <- make.predictorMatrix(prepped_df)
predMat[,c(1:2,22:24)] <- 0
meth <- make.method(prepped_df)

# impute
# 20% missing, then use 20 imputed datasets
imputed_df <- mice(prepped_df,m=5,maxit=50,meth=meth,seed=500, predictorMatrix=predMat)
summary(imputed_df)
completed_df <- complete(imputed_df,1)

######################################################################
###### run mixed graphical model network for mixed data types ########
###################################################################### 

all_network_list <- list()

# first estimate the network using
for (wave in c(1,3,5)) {
  timepoint <- completed_df$time == wave
  data_subset <- completed_df[timepoint, 3:ncol(completed_df)]
  all_network <- estimateNetwork(data_subset, default="mgm")
  all_network_list[[as.character(wave)]] <- all_network
}

#class_values <- c("Decreasing", "Stable low", "Persistent", "Increasing")

'for (wave in c(1,3,5,7)) {
  for (class_value in class_values) {
    timepoint <- completed_df$time == wave
    class_condition <- symp_gen_class$class == class_value
    data_subset <- completed_df[timepoint #& class_condition
                                  , 3:ncol(completed_df)]
    all_network <- estimateNetwork(data_subset, default = "mgm")
    all_network_list[[paste0("Wave", wave
                             , "_Class", class_value
                             )]] <- all_network
  }
}'

# network colours: can be "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
theme <- 'colorblind'
# node colours: can be "rainbow","colorblind", "pastel","gray","R","ggplot2"
palette <- 'colorblind'


groups=list(1:13,14:19,20:22)
group_names <- c("SMFQ symptoms", "Environmental factors", "Genetic factors")
names(groups) <- group_names

# Loop through the results and generate plot for estimated networks
for (wave in names(all_network_list)) {
  results <- all_network_list[[wave]]
  qgraph(results$graph, 
         layout='spring',
          theme=theme, groups=groups,
         palette = palette, nodeNames=c(unlist(labels),c(unlist(env_labels),'mood','psychotic','neurodevelopmental')),
         legend=TRUE, legend.mode='style1', legend.cex=0.4)
  title(paste(wave), adj=0.1, line=-0.8)
  #qgraph::centralityPlot(results, include = c("Strength","Closeness","Betweenness"), 
             #labels=c(unlist(labels),c(unlist(env_labels), 'mood','psychotic','neurodevelopmental')),
              #orderBy = 'default', theme_bw = TRUE, scale = 'z-scores')
  }

## overlay centrality plot
centralityPlot(list(Wave1 = all_network_list[[1]], Wave3 = all_network_list[[2]], 
                    Wave5 = all_network_list[[3]]), 
               theme_bw=FALSE, scale = "z-scores", 
               include = c("Strength","Closeness","Betweenness"), 
               labels = nodeNames)


## 2: Run bootnet resampled network
all_boot_result_list <- list()

## for each wave run bootnet and store results in list 
for (wave in names(all_network_list)) {
  results <- all_network_list[[wave]]
  boot1<-bootnet(results, nBoots=50, default="mgm", nCores=1)
  all_boot_result_list[[as.character(wave)]] <- boot1
}

## plot bootstrapped network results
for (wave in names(all_boot_result_list)) {
  results <- all_boot_result_list[[wave]]
  plot(results$sample, label = c(unlist(labels),c(unlist(env_labels),'mood','psychotic','neurodevelopmental')), theme=theme)
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
grid.arrange(grobs = edge_plot_list, nrow=1, ncol=3)

## Plot significant differences (alpha=0.05) of nodestrength
# grey = non sig, black = sig, white = value of node strength 
# Create a list of strength difference plots and arrange in 1xn grid
strength_plot_list <- lapply(all_boot_result_list, function(boot_result) {
  plot(boot_result, "strength", plot = "difference", order = "sample",
       labels=TRUE, legend=TRUE)})
grid.arrange(grobs = strength_plot_list, nrow=1, ncol=3)

########################################
########## symptom heatmaps ###########

nodeNames = c(unlist(labels),c(unlist(env_labels),'mood','psychotic','neurodevelopmental'))

result1 <- all_network_list[["1"]][["graph"]]
result3 <- all_network_list[["3"]][["graph"]]
result5 <- all_network_list[["5"]][["graph"]]

results <- list(result1, result3, result5)

for (i in seq_along(results)) {
  rownames(results[[i]]) <- nodeNames
  colnames(results[[i]]) <- nodeNames
}

for (i in seq_along(results)) {
  heatmap(results[[i]], 
          symm = TRUE,
          col = viridis::plasma(100),
          Rowv = NA,
          main = paste("Wave", i))
}


#####################################

## Network comparison test
wave1 <- completed_df %>% filter(time == 1)
x <- estimateNetwork(wave1[,3:24], default='mgm')
wave3 <- completed_df %>% filter(time == 3)
y <- estimateNetwork(wave3[,3:24], default='mgm')
wave5 <- completed_df %>% filter(time == 5)
z <- estimateNetwork(wave3[,3:24], default='mgm')

NCT(x, y)
NCT(y, z)

#####################################
