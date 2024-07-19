#############################################################################
############################# MCS network temp ##############################
#############################################################################
library(foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library(qgraph)
library(mice)
library(ggmice)
library(VIM)

############################# LOAD AND QC ######################################
# read in data and label cols
setwd('/Volumes/igmm/GenScotDepression/users/poppy/mcs/symptom_data')
mcs_qcd <- read.table('mcs_sdq_sym_long.txt', check.names = FALSE)

labels <- c("malaise", "worries", "unhappy", "anxiety","fears", 
            "solitary","friends","liked","bullied", "adultoriented")

emo.labels <- c("malaise", "worries", "unhappy", "anxiety","fears")

colnames(mcs_qcd) <- c('id', 'sex', 'time', unlist(labels))

# recode variables [-1,1] encoding
# 'friends' and 'liked' are reverse coded
mcs_qcd <- mcs_qcd %>%
  mutate(across(c(4:9, 12:13), ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1,
    TRUE ~ .
  ))) %>%
  mutate(across(10:11, ~case_when( # liked and friends reverse scored
    . == 1 ~ -1,
    . == 0 ~ 1,
    TRUE ~ .
  )))

# select emotional scale symptoms only
mcs.emo <- mcs_qcd[1:8]

################################ IMPUTATION ####################################
# prepare the df and merge with original data
template <- expand.grid(id = unique(mcs.emo$id), time = 5:7) # try 4:7 and 5:7?
prepped_df <- merge(template, mcs.emo, by = c("id", "time"), all.x = TRUE)
# sort by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(`malaise`:`fears`, ~replace(., is.na(.), NA)))
# fill sex (have checked that individuals don't have conflicting sex)
prepped_df <- prepped_df %>%
  group_by(id) %>%
  fill(sex, .direction = "down") %>% fill(sex, .direction = "up")
# check prepped df
print(prepped_df)

### inspect missing data 
plot_pattern(prepped_df, square=FALSE) + 
  theme(plot.title = element_text(hjust = 0.5)) + ggtitle("MCS")
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
## A lot of missing sex data in MCS? N=8776

### make predictor matrix (but don't impute time or sex)
predMat <- make.predictorMatrix(prepped_df)
predMat[,1] <- 0
meth <- make.method(prepped_df)
meth[c("sex")] <- ""

#### impute (if error, check var labels)
imputed_mcs <- mice(prepped_df, m=33 ,maxit=20 ,meth=meth, seed=123, predictorMatrix=predMat)
saveRDS(imputed_mcs, 'imputed_mcs_33.rds')

# extract eand combine pooled datasets into long format
imputed_list <- lapply(1:33, function(i) complete(imputed_mcs, action = i))
all_imputed <- bind_rows(imputed_list)

# sanity check: count the number of duplicate rows based on 'id' and 'time' (14,958*4 = 59832)
nrow(all_imputed %>%
  group_by(id, time, sex) %>%
  filter(n() > 1) %>%
  summarise(count = n()) %>%
  ungroup())

# group by 'id', 'time', 'sex' and the symptom columns, then count occurrences
combination_counts <- all_imputed %>%
  group_by(id, time, sex, malaise, worries, unhappy, anxiety, fears) %>%
  summarise(count = n(), .groups = 'drop')

# for each grouping of 'id', 'time', 'sex', select the combination with the highest count
mcs_mode_imputed <- combination_counts %>%
  group_by(id, time, sex) %>%
  slice_max(count, n = 1, with_ties = FALSE) %>%
  ungroup() %>% dplyr::select(-count)

# final dataset that pools the mode symptom combination from multiple imputation
head(mcs_mode_imputed) ## nrow = 59832

################################################################################
### sex stratification
boys <- mcs_mode_imputed %>% filter(sex==0)
girls <- mcs_mode_imputed %>% filter(sex==1)

## [0,1] encoding
imputed_mcs_recode <- mcs_mode_imputed %>%
  mutate(across(4:8, ~case_when(
    . == 1 ~ 1,
    . == -1 ~ 0,
    TRUE ~ .
  )))

########################## MULTIGROUP ISING MODEL ##############################

# symptom vars to include in network
vars <- emo.labels

# remove age 7 (t4)
no.t4.og <- mcs.emo[mcs.emo$time != 4, ] # without imputation (complete cases)
#no.t4.imp <- mcs_mode_imputed[mcs_mode_imputed$time != 4, ] # after imputation
boys.no.t4 <- boys[boys$time != 4,]
girls.no.t4 <- girls[girls$time != 4,]

# Form saturated model and run [all params free]
model1 <- Ising(boys, vars = vars, groups = "time") %>% runmodel
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

# Compare all models, increasing constraints (RMSEA <0.05 good, minimise AIC and BIC)
psychonetrics::compare(
  `1. all parameters free (dense)` = model1, # saturated
  `2. all parameters free (sparse)` = model1b,
  `3. equal networks (dense)` = model2, # equal network structure 
  `4. equal networks (sparse)` = model2b,
  `5. equal networks and thresholds (dense)` = model3, # structure and external fields
  `6. equal networks and thresholds (sparse)` = model3b,
  `7. all parameters equal (dense)` = model4, # structure, external and temperature 
  `8. all parameters equal (sparse)` = model4b
) %>% arrange(BIC) 

best.model <- model2
  
# extract and plot network
# "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
network_sdq <- getmatrix(best.model, "omega")[[1]]
graph_sdq <- qgraph(network_sdq, layout = 'spring', labels = vars, 
                     theme = 'colorblind', label.prop=0.99, node.width=1.4, 
                     label.norm='000000')

# centrality plot
all_network_sdq <- getmatrix(best.model, "omega")
centralityPlot(list(Wave1 = all_network_sdq[[1]]),
theme_bw=FALSE, scale = "z-scores", 
include = c("Strength","Closeness","Betweenness"), 
labels = vars, orderBy = 'Strength')

# extract temperature and calculate 95%CIs from SE of beta
temp_sdq <-  as.numeric(lapply(getmatrix(best.model, "beta"), 'mean'))
betas_est <- best.model@parameters$est[best.model@parameters$matrix == "beta"]
betas_se <- best.model@parameters$se[best.model@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_sdq + (z*betas_se)
lowerCI <- temp_sdq - (z*betas_se)

# extract external fields
fields_sdq <- lapply(getmatrix(best.model, 'tau'), 'mean')

# plot heatmap
rownames(network_sdq) <- emo.labels
colnames(network_sdq) <- emo.labels

heatmap(network_sdq, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA,
        main = paste("SDQ items"))

# gather temp data
ages <- c(11,14,17)
mcs_temp_data_boys <- data.frame(Age = ages,
                        Temperature = 1/temp_sdq,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)

mcs_temp_data_girls$sex <- "Females"
mcs_temp_data_boys$sex <- "Males"
mcs_temp_data <- rbind(mcs_temp_data_boys, mcs_temp_data_girls)

# test whether rate of decrease is sig different
# look at interaction term age:sex
model <- lmer(Temperature ~ Age * sex + (1 | Age), data = mcs_temp_data)
summary(model)

# plot overlay temperature change and confidence intervals
ggplot(mcs_temp_data, aes(x = Age, y = Temperature, color=sex)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.6, 1.0) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),  
        plot.title = element_text(hjust = 0.5, size = 12), legend.position='none') +  # center and adjust size of title
  scale_x_continuous(breaks = ages) +
  scale_color_manual(values = c("Females" = "darkgrey", "Males" = "black"))

# panel plot
ggplot(temp_data, aes(x = Age, y = Temperature, group = sex)) +
  geom_line(color = "black") +                    
  geom_point(color = "black") +     
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") + 
  facet_wrap(~ sex) +                        
  theme_bw() +    
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = 'NULL',
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text = element_text(size = 16)) +
  labs(x = "Age", y = "Temperature")

####################### BOOTSTRAPPING TEMPERATURE #########################

# read in results from eddie/datastore
setwd('/Volumes/igmm/GenScotDepression/users/poppy/mcs')
boot.mcs <- readRDS('bootstrap_results_mcs.rds')
print(boot.mcs)

std_errors <- c(NA,apply(boot.mcs$t, 2, sd))
z = qnorm(0.975)
upper_boot_CI <- temp_sdq + (z*std_errors)
lower_boot_CI <- temp_sdq - (z*std_errors)

ages <- c(11,14,17)
boot_ci_data <- data.frame(Age = ages,
                           Temperature = 1/temp_sdq,
                           UpperCI = 1/upperCI,
                           LowerCI = 1/lowerCI,
                           UpperBootCI = 1/upper_boot_CI,
                           LowerBootCI = 1/lower_boot_CI)


boot_ci_data <- boot_ci_data %>%
  mutate(CI_Type = 'Analytical 95% CIs',
         BootCI_Type = 'Bootstrapped 95% CIs')

ggplot(boot_ci_data, aes(x = Age, y = Temperature)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI, color = CI_Type), width = 0.1, show.legend = TRUE) +
  geom_errorbar(aes(ymin = LowerBootCI, ymax = UpperBootCI, color = BootCI_Type), width = 0.1, show.legend = TRUE) +
  ylim(0.58, 1.0) +
  labs(x = "Age (years)", y = "Temperature") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = 'NULL',
    axis.ticks.length = unit(14, "points")) +
  scale_x_continuous(breaks = ages) +
  scale_color_manual(values = c('Analytical 95% CIs' = 'purple', 'Bootstrapped 95% CIs' = 'red'), name=NULL) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "solid"))))


##################### SENSITIVITY/SUPPLEMENTAL ANALYSES ########################

# histograms of overall depression score 
mcs_qcd$total <- rowSums(mcs_qcd[,3:12])
breaks <- seq(min(mcs_qcd$total, na.rm = TRUE), max(mcs_qcd$total, na.rm = TRUE), length.out = 10)
hist(mcs_qcd$total[mcs_qcd$time==5], main = 'Age 11', xlab = 'Overall depression', breaks=breaks)
mtext('(c)', 3, at = -46, padj = -4)
par(mar = c(6, 3, 6, 3))
hist(mcs_qcd$total[mcs_qcd$time==6], main = 'Age 14', xlab = 'Overall depression', breaks=breaks)
par(mar = c(6, 2, 6, 4))
hist(mcs_qcd$total[mcs_qcd$time==7], main = 'Age 17', xlab = 'Overall depression', breaks=breaks)
dev.off()

## summary score distribution
# add sum score col
mcs_mode_imputed <- mcs_mode_imputed %>%
  dplyr::mutate(total = rowSums(dplyr::select(., malaise:fears)))

mcs_variance <- mcs_mode_imputed %>%
  group_by(time) %>%
  summarise(variance = var(total, na.rm = TRUE))

print(mcs_variance)

# histogram of sumscore (can plot for [-1,1] and [0,1] encoding)
ggplot(mcs_mode_imputed, aes(x = total)) +
  geom_histogram(binwidth = 2, fill = "skyblue", color = "black") +
  facet_wrap(~time, ncol = 3) +
  labs(title = "Distribution of Sumscore by Time",
       x = "Sumscore",
       y = "Frequency")

####### variance of original sample and complete cases [for Table 1]

# labelled t4-6 but have checked ages 11,14,17 
mcs.sumscores <- read.csv('/Volumes/igmm/GenScotDepression/users/poppy/mcs/4_sdq_waves.csv') %>%
  dplyr::select(c(1,4:6))

colnames(mcs.sumscores) <- c('id','t1','t2','t3')

# complete cases variance
complete_cases <- mcs.sumscores[complete.cases(mcs.sumscores), ]
summary_stats <- complete_cases %>%
  summarise(across(t1:t3, list(mean = mean, SD = sd), na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = c("Variable", ".value"), names_sep = "_")
print(summary_stats)

mcs.long <- complete_cases %>% 
    pivot_longer(
        cols = `t1`:`t3`, 
       names_to = "time",
       values_to = "value"
     )

ggplot(mcs.long, aes(x = value)) +
  geom_histogram(binwidth = 2, fill = "skyblue", color = "black") +
  facet_wrap(~time, ncol = 3) +
  labs(title = "Distribution of Sumscore by Time",
       x = "Sumscore",
       y = "Frequency")
#############################################

