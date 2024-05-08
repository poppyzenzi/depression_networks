#################################### 
######## ALSPAC network temp #######
#################################### 
library(foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library (qgraph)
library(mice)
library(ggmice)
library(gplots)
library(VIM)
library(boot)
library(purrr)
################## DATA AND QC ####################

## read in the data, add predictor variables and labels
setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
smfq_symptoms <- read.table('smfq_symptoms_qcd.txt', check.names = FALSE)
env <- read.table('alspac_envi_vars.txt', check.names=FALSE)
smfq_qcd <- merge(smfq_symptoms, env, by='id')
smfq_qcd <- smfq_qcd %>% rename(`sex` = kz021,
                                `maternal_depression` = r2021,
                                `bullying` = f8fp470,
                                `child_trauma` = AT5_n,
                                `sleep` = FJCI250,
                                `income` = h470,
)

env_labels = list(names(smfq_qcd[,16:ncol(smfq_qcd)]))
labels <- c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
            "tearful", "distracted", "self_loathing", "guilty", "isolated", 
            "unloved", "inadequate", "incompetent")
colnames(smfq_qcd) <- c('id', 'time', unlist(labels), unlist(env_labels))

# recode variables so that each symptoms is binary with +1 and -1 
smfq_qcd <- smfq_qcd %>%
  mutate(across(3:15, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1,
    TRUE ~ .
  )))

########## IMPUTATION ##############

### Prepare the df and merge with og df
template <- expand.grid(id = unique(smfq_qcd$id), time = 1:7)
prepped_df <- merge(template, smfq_qcd, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(unhappy:incompetent, ~replace(., is.na(.), NA)))
# fill sex (have checked that individuals don't have conflicitng sex)
prepped_df <- prepped_df %>%
  group_by(id) %>%
  fill(sex, .direction = "down") %>% fill(sex, .direction = "up")
# reorder cols
prepped_df <- prepped_df[, c("id", "time", unlist(labels), unlist(env_labels))]
# Check completed df
print(prepped_df)

### inspect 
plot_pattern(prepped_df, square = TRUE, rotate = TRUE) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ALSPAC")
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

### make predictor matrix
# we only want to impute the symptoms vars for now
predMat <- make.predictorMatrix(prepped_df)
predMat[,c(1:2,16:21)] <- 0
meth <- make.method(prepped_df)
meth[c("sex","maternal_depression","bullying","child_trauma","sleep","income")] <- ""

#### impute
imputed_smfq <- mice(prepped_df,m=5,maxit=40,meth=meth,seed=500, predictorMatrix=predMat)
summary(imputed_smfq)
imputed_df_smfq <- complete(imputed_smfq,1)

write.table(imputed_df_smfq, 'alspac_imputed_ising.txt', col.names=TRUE)

########## options to stratify the data ###########

# reappend the env vars 
# stratify by sex
girls <- subset(imputed_df_smfq, sex==1)
boys <- subset(imputed_df_smfq, sex==0)
mat_dep <- subset(imputed_df_smfq, maternal_depression==1)
no_mat_dep <- subset(imputed_df_smfq, maternal_depression==0)
######### fit psychonetrics Ising model ##########

#fit an Ising model with increasing constraints representing their hypotheses to this longitudinal assessment 
# We investigate the impact on the fit of the model of 
# 1) constraining edges between nodes to be equal across time points
# 2) constraining the external fields to be equal across time points 
# 3) constraining the temperature (the entropy of the system) to be equal across time points. 
# 4) Additionally, we tested whether a dense network (all nodes are connected) or 
#  a sparse network (at least some edges are absent) fits the data best

# Variables to use:
vars <- names(imputed_df_smfq)[3:15]

# remove age 22
no.t7 <- imputed_df_smfq[imputed_df_smfq$time != 7, ] %>% as.data.frame()
boys.no.t7 <- no.t7 %>% filter(sex==0)
girls.no.t7 <- no.t7 %>% filter(sex==1)

# Form saturated model and run [all params free]
model1 <- Ising(girls.no.t7, vars = vars, groups = "time")
## here add estimator='FIML' to the Ising model for imputation
model1 <- model1 %>% runmodel
# Prune-stepup to find a sparse model:
model1b <- model1 %>% prune(alpha = 0.05) %>%  stepup(alpha = 0.05)
# Equal networks (omega = network structure, edges btw nodes equal across time):
model2 <- model1 %>% groupequal("omega") %>% runmodel
# Prune-stepup to find a sparse model:
model2b <- model2 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)
# Equal thresholds (tau = threshold/intercept structure, omega still equal plus external fields equal):
model3 <- model2 %>% groupequal("tau") %>% runmodel
# Prune-stepup to find a sparse model:
model3b <- model3 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)
# Equal beta (beta = inverse temperature, omega and tau also constrained):
model4 <- model3 %>% groupequal("beta") %>% runmodel
# Prune-stepup to find a sparse model:
model4b <- model4 %>% prune(alpha = 0.05) %>% stepup(mi = "mi_equal", alpha = 0.05)

# Compare all models, increasing constraints:
# RMSEA <0.05 good, minimise AIC and BIC
comparison <- psychonetrics::compare(
    `1. all parameters free (dense)` = model1, # saturated
    `2. all parameters free (sparse)` = model1b,
    `3. equal networks (dense)` = model2, # equal network structure 
    `4. equal networks (sparse)` = model2b,
    `5. equal networks and thresholds (dense)` = model3, # structure and external fields
    `6. equal networks and thresholds (sparse)` = model3b,
    `7. all parameters equal (dense)` = model4, # structure, external and temperature 
    `8. all parameters equal (sparse)` = model4b
  ) %>% arrange(BIC) 

print(comparison)
best.model <- model2

#extract and plot network
all_network_smfq <- getmatrix(best.model, "omega")
network_smfq <- getmatrix(best.model, "omega")[[1]]
graph_smfq <- qgraph(network_smfq, layout = 'spring', labels = vars, theme = 'colorblind')

## overlay centrality plot
centralityPlot(list(Wave1 = all_network_smfq[[1]]),
               theme_bw=FALSE, scale = "z-scores", 
               include = c("Strength","Closeness","Betweenness"), 
               labels = vars, orderBy = 'Strength')

#extract temperature
temp_smfq <- as.numeric(lapply(getmatrix(best.model, "beta"), 'mean'))
# calculate 95%CIs from standard errors of the beta parameter
betas_est <- best.model@parameters$est[best.model@parameters$matrix == "beta"]
betas_se <- best.model@parameters$se[best.model@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_smfq + (z*betas_se)
lowerCI <- temp_smfq - (z*betas_se)

## plot temperature change
ages <- c(11,13,14,17,18,19)
income <- 1:5
temp_data_girls <- data.frame(Age = ages,
                        Temperature = 1/temp_smfq,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)

temp_data_girls$Sex <- "Females"
temp_data_boys$Sex <- "Males"
temp_data <- rbind(temp_data_boys, temp_data_girls)
temp_data_matdep$matdep <- "Mat dep"
temp_data_nomatdep$matdep <- "No mat dep"
temp_data <- rbind(temp_data_matdep, temp_data_nomatdep)

ggplot(temp_data, aes(x = Age, y = Temperature, color=Sex)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.8, 1.0) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), 
        plot.title = element_text(hjust = 0.5, size = 12)) + 
  scale_x_continuous(breaks = ages) + 
  scale_color_manual(values = c("Females" = "orange", "Males" = "blue"))


#extract external fields
fields_smfq <- lapply(getmatrix(model2, 'tau'), 'mean')

######################################################


######################################################
################### bootstrap temp ###################

## 1. data, 2. statistic (temperature), 3. run bootstrapping with N boots, 4. get bootstrapped CIs

# 1. check the data (longitudinal binary symptoms +1 -1 encoding)
head(imputed_df_smfq)
vars <- names(imputed_df_smfq[,3:15])

# 2. statistic function to extract a vector of temperatures
temperature_statistic <- function(data, indices) {
    boot_data <- data[indices, ]
    # run Ising network model on boot data
    model <- Ising(boot_data, vars=vars, groups = 'time') %>% groupequal("omega") %>% runmodel
    # extract temperature estimates from model (1/beta)
    temp_est <- 1/model@parameters$est[model@parameters$matrix=='beta']
    # omit t1 as fixed to 1 and stops boot.ci running
    temp_est <- temp_est[-1]
    return(temp_est)
} 

# 3. perform bootstrap resampling
bootstrap_results <- boot(data = imputed_df_smfq, statistic = temperature_statistic, R = 1000)
bootstrap_results2 <- boot(data = imputed_df_smfq, statistic = temperature_statistic, R = 3)

# 4a. get CIs
boot.ci(bootstrap_results, type='bca')
boot.ci(bootstrap_results2, type='bca')
# get the output bootstrap statistics
original <- bootstrap_results$t0
temperatures <- bootstrap_results$t
bias <- mean(bootstrap_results$t)-bootstrap_results$t0
SE <- sd(bootstrap_results$t)

# 4b. get CIs manually

setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/projects/networks/depression_networks')
## run the above steps 1-3 with 100 and 1000 bootstraps on eddie
boots100 <- readRDS('from_eddie/bootstrap_R100.rds')
boots1k <- readRDS('from_eddie/bootstrap_R1000.rds') # one thousand boots

# number of estimates (groups)
num_columns <- ncol(boots1k$t)

# vectors to store the quantiles
quantile_975 <- numeric(num_columns)
quantile_025 <- numeric(num_columns)

# calculate quantiles for each estimate
for (i in 1:num_columns) {
  quantile_975[i] <- quantile(boots1k$t[, i], probs = 0.975)
  quantile_025[i] <- quantile(boots1k$t[, i], probs = 0.025)
}

cat("97.5th percentile:", quantile_975, "\n")
cat("2.5th percentile:", quantile_025, "\n")

# histogram of bootstrap temperatures
hist(boots1k$t, main = "Bootstrap Distribution of Temperature", xlab = "Temperature", 
     ylab = "Frequency", col = "skyblue", border = "white")

## plot temperature change
ages <- c(11,13,14,17,18,19,22)
# plot estimates
boot_temp_data <- data.frame(Age = ages,
                        Temperature = c(1, boots1k$t0),
                        bootUpperCI = c(1,quantile_975),
                        bootLowerCI = c(1,quantile_025))

# Plot estimates and confidence intervals using ggplot2
ggplot(boot_temp_data, aes(x = Age, y = Temperature)) +
  geom_point(color = "blue", shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = bootLowerCI, ymax = bootUpperCI), width = 0.2) +
  #ylim(0.80, 1.1) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),  
        axis.text.y = element_text(size = 11),  
        plot.title = element_text(hjust = 0.5, size = 12)) +  # Center and adjust size of title
  scale_x_continuous(breaks = ages)

## dataset being bootstrapped at X-sectional time points treating each sample independently
## Doesn't account for correlations between measurements from same individuals over time?
## Numbers too small?
## Binary encoding can't capture variation within the data?
## Check step 2 all ok?


# bootstrap IDs
###########################################

# Imputed data: remove sex, convert to wide, select vars 
imputed_df_smfq <- imputed_df_smfq[, !(names(imputed_df_smfq) == "sex")] 
wide_data <- reshape(imputed_df_smfq, idvar = "id", timevar = "time", direction = "wide") 

# Function to take wide dataset and bootstrap on individuals 
bootstrap_function <- function(data, indices) {
  boot_data <- data[indices, ]
  boot_data$id <- seq(1, nrow(boot_data))
  # Convert wide dataset to long format
  boot_data_long <- reshape(boot_data, 
                            varying=list(unhappy = seq(2, 92, by = 13),
                                         anhedonia = seq(3, 93, by = 13),
                                         apathetic = seq(4, 94, by = 13),
                                         restless = seq(5, 95, by = 13),
                                         worthless = seq(6, 96, by = 13),
                                         tearful = seq(7, 97, by = 13),
                                         distracted = seq(8, 98, by = 13),
                                         self_loathing = seq(9, 99, by = 13),
                                         guilty = seq(10, 100, by = 13),
                                         isolated = seq(11, 101, by = 13),
                                         unloved = seq(12, 102, by = 13),
                                         inadequate = seq(13, 103, by = 13),
                                         incompetent = seq(14, 104, by = 13)),
                            v.names=c("unhappy","anhedonia","apathetic","restless","worthless","tearful","distracted","self_loathing" ,"guilty",
                                      "isolated","unloved","inadequate","incompetent"),
                            direction="long",
                            time=1:7,
                            timevar="time")
  # Select variables for network analysis
  vars <- names(imputed_df_smfq)[3:15]
  # Perform network analysis on boot data in long format
  model <- Ising(boot_data_long, vars = vars) %>% groupequal("omega") %>% runmodel
  # Extract temperature estimates from model (1/beta)
  temp_est <- 1 / model@parameters$est[model@parameters$matrix == 'beta']
  # Omit t1 as it's fixed to 1 and stops boot.ci running
  temp_est <- temp_est[-1]
  return(temp_est)
}

# Number of bootstrap samples
num_bootstraps <- 5

# Run bootstrap
bootstrap_results <- boot(data = wide_data, statistic = bootstrap_function, R = num_bootstraps)



################################################
################### plotting ################### 


#plot heatmaps of symptom covar matrix
rownames(network_smfq) <- labels[1:13]
colnames(network_smfq) <- labels[1:13]

heatmap(network_smfq, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA,
        main = paste("SMFQ items"))

heatmap.2(network_smfq, 
          symm = TRUE,
          col = viridis::plasma(100),Rowv = NA,
          trace = "none", density.info = "none")

par(mar = rep(0,4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
                   "tearful", "distracted", "self-loathing", "guilty", "isolated", 
                   "unloved", "inadequate", "incompetent"), 
       title = 'SMFQ items', pch = 19,
       cex = 0.8, bty = 'n')

# a label at 3rd margin (top) at 0.1 along and right justification (2)
mtext('(a)', 3, at = .01, padj = 2)

# palette (for nodes): "rainbow"(default), "colorblind","pastel","gray","R"and"ggplot2"
# theme (for edges): "classic", "colorblind", "gray", "Hollywood", "Borkulo", "gimme", "TeamFortress", "Reddit", "Leuven" or "Fried".

qgraph(network_smfq, layout = 'spring', 
       groups = list(1:13), 
       posCol='slateblue', negCol='red', palette='colorblind',
       legend = FALSE,
       labels = c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
                  "tearful", "distracted", "self-loathing", "guilty", "isolated", 
                  "unloved", "inadequate", "incompetent"), vsize = 10)

par(mar = rep(2,4), cex.main = 0.8)

## plot temperature change
ages <- c(11,13,14,17,18,19,22)
# plot estimates
temp_data <- data.frame(Age = ages,
                   Temperature = 1/temp_smfq,
                   UpperCI = 1/upperCI,
                   LowerCI = 1/lowerCI)

# Plot estimates and confidence intervals using ggplot2
ggplot(temp_data, aes(x = Age, y = Temperature)) +
  geom_point(color = "blue", shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.80, 1.0) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),  # Adjust size of x-axis tick labels
        axis.text.y = element_text(size = 11),  # Adjust size of y-axis tick labels
        plot.title = element_text(hjust = 0.5, size = 12)) +  # Center and adjust size of title
  scale_x_continuous(breaks = ages)


# histograms of overall depression score 
imputed_df_smfq_smfq$total <- rowSums(imputed_df_smfq_smfq[,3:15])

hist(imputed_df_smfq_smfq$total[imputed_df_smfq_smfq$time==1], main = 'Age 11', xlab = 'Overall depression')
mtext('(c)', 3, at = -46, padj = -4)

# Set up the layout for multiple plots
par(mfrow=c(2, 4))

# Plot each histogram
hist(imputed_df_smfq$total[imputed_df_smfq$time==1], main = 'Age 11', xlab = 'Overall depression')
hist(imputed_df_smfq$total[imputed_df_smfq$time==2], main = 'Age 13', xlab = 'Overall depression',  yaxt = 'n', ylab = '')
hist(imputed_df_smfq$total[imputed_df_smfq$time==3], main = 'Age 14', xlab = 'Overall depression',  yaxt = 'n', ylab = '')
hist(imputed_df_smfq$total[imputed_df_smfq$time==4], main = 'Age 17', xlab = 'Overall depression',  yaxt = 'n', ylab = '')
hist(imputed_df_smfq$total[imputed_df_smfq$time==5], main = 'Age 18', xlab = 'Overall depression',  yaxt = 'n', ylab = '')
hist(imputed_df_smfq$total[imputed_df_smfq$time==6], main = 'Age 19', xlab = 'Overall depression',  yaxt = 'n', ylab = '')
hist(imputed_df_smfq$total[imputed_df_smfq$time==7], main = 'Age 22', xlab = 'Overall depression',  yaxt = 'n', ylab = '')


plot_images <- list()

# Loop through each time point and create a histogram plot
for (i in 1:7) {
  # Create histogram plot
  hist_plot <- hist(imputed_df_smfq$total[imputed_df_smfq$time==i], 
                    main = paste("Age", 10 + i), 
                    xlab = "Overall depression", 
                    yaxt = "n", 
                    ylab = "")
  
  # Convert the plot to an image
  plot_image <- image_graph(width = 800, height = 600)
  plot(hist_plot)
  dev.off()
  
  # Store the image in the list
  plot_images[[i]] <- plot_image
}

# Combine the images into a GIF
gif_file <- image_animate(image_join(plot_images), fps = 1)

# Save the GIF file
gif_file_path <- "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/projects/networks/depression_networks/figures/histograms.gif"
image_write(gif_file, gif_file_path)

hist(smfq_qcd$total[smfq_qcd$time==4], main = 'Age 17', xlab = 'Overall depression',  yaxt = 'n', ylab = '')

dev.off()


### Compare networks of data sets using NCT ###
## Networks can be compared by either (1) feeding the data directly into NCT (whereby 
## you need to specify arguments such as "gamma" and "binary.data") or (2) by using 
## estimateNetwork() (bootnet package) and feeding that output into NCT. For the latter 
## option, we refer to the help file of  estimateNetwork() for its usage. Below, both 
## options are illustrated. We recommend using estimateNetwork(), since this function 
## has implemented many network estimation methods.

## gamma = 0 (in estimateNetwork this hyperparameter is called "tuning"; to illustrate 
# how to specify a different value than the default)
## iterations (it) set to 10 to save time
## Note: Low number of iterations can give unreliable results; should be 1000 at least

## Testing whether there are differences in the three aspects that are validated 
# (network invariance, global strength, edge weight)
## 2 edges are tested here: between variable 1 and 2, and between 3 and 6 (can be 
# "list(c(2,1),c(6,3))" as well)

## (1) Feeding data directly into NCT, plot results of global strength invariance test
set.seed(123)
NCT_a <- NCT(subset(smfq_qcd, time==1)[,3:15], subset(smfq_qcd, time==2)[,3:15], gamma=0, it=50, binary.data = TRUE, 
             test.edges=TRUE, edges=list(c(1,2),c(3,6)))
summary(NCT_a)
plot(NCT_a, what="strength")

NCT_b <- NCT(subset(smfq_qcd, time==2)[,3:15], subset(smfq_qcd, time==3)[,3:15], gamma=0, it=50, binary.data = TRUE, 
             test.edges=TRUE, edges=list(c(1,2),c(3,6)))
summary(NCT_b)
plot(NCT_b, what="strength")

NCT_c <- NCT(subset(smfq_qcd, time==3)[,3:15], subset(smfq_qcd, time==4)[,3:15], gamma=0, it=50, binary.data = TRUE, 
             test.edges=TRUE, edges=list(c(1,2),c(3,6)))
summary(NCT_c)
plot(NCT_c, what="strength")

#####


# genetic quintile testing


# genetically straitfy by quintiles
quintiles <- quantile(smfq_qcd$MDDPRS, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1))
labels <- c("very low", "low", "middle", "high", "very high")
smfq_qcd$quintile <- cut(smfq_qcd$MDDPRS, breaks = quintiles, labels = labels, include.lowest = TRUE)

# get distribution and frequency
hist(smfq_qcd$mood)
table(smfq_qcd$quintile)

# get mean and SD of PRS in each group 
risk_group_types <- unique(smfq_qcd$quintile)  # Get quintiles

very_high <- subset(smfq_qcd, quintile=='very high')
very_low <- subset(smfq_qcd, quintile=='very low')
low <- subset(smfq_qcd, quintile=='low')
high <- subset(smfq_qcd, quintile=='high')
med <- subset(smfq_qcd, quintile=='middle')
high_quart <- subset(smfq_qcd, HighQuartile==1) # nunique = 1552

# dataframe for summary statistics
summary_statistics <- data.frame(Category = character(),
                                 Mean = numeric(),
                                 SD = numeric(),
                                 Count = integer(),  # Change to integer type
                                 stringsAsFactors = FALSE)
# Loop through each risk group type
for (risk_group_type in risk_group_types) {
  mean_value <- mean(smfq_qcd$MDDPRS[smfq_qcd$risk_group == risk_group_type], na.rm = TRUE)
  sd_value <- sd(smfq_qcd$MDDPRS[smfq_qcd$risk_group == risk_group_type], na.rm = TRUE)
  count_value <- count[risk_group_type]
  summary_statistics <- rbind(summary_statistics, 
                              data.frame(Category = risk_group_type, 
                                         Mean = mean_value, 
                                         SD = sd_value, 
                                         Count = count_value))
}

summary_statistics

## highest quartile method
quartiles <- quantile(smfq_qcd$MDDPRS, probs = c(0.8))  # 75th percentile (highest quartile)
smfq_qcd$HighQuartile <- ifelse(smfq_qcd$MDDPRS >= quartiles, 1, 0)
head(smfq_qcd)