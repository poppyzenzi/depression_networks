#############################################################################
############################# ABCD network temp ##############################
#############################################################################
library (foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library (qgraph)
library(VIM)
library(mice)
library(ggmice)

############################# LOAD AND QC ######################################
# data and labels (symptoms at t1:t8)
setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")
abcd_qcd <- read.table('abcd_bpm_sym_long.txt', check.names = FALSE)
labels <- c("worthless","anxious","guilty","self_conscious","unhappy","worry")
colnames(abcd_qcd) <- c('id', 'time', unlist(labels))

# demographics
demog <- read.table('/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_envi_vars.txt')
env.vars <- demog %>% filter(eventname==0) %>% dplyr::select(src_subject_id, demo_sex_v2, fam_history_q6d_depression, demo_comb_income_v2)
names(env.vars) <- c('id','sex','matdep','income')
abcd_qcd <- merge(abcd_qcd, env.vars, by='id', all.x=TRUE)

# puberty (at t0,2,4,6,8, t0 is lost as no symptoms)
puberty <- read.table('../predictors/puberty.vars.txt') %>% dplyr::select(src_subject_id, eventname, pds_y_ss)
names(puberty) <- c('id','time','pubertal_stage')
abcd_qcd <- merge(abcd_qcd, puberty, by=c('id','time'), all.x=TRUE)

# recode variables [-1,1]
abcd_qcd <- abcd_qcd %>%
  mutate(across(3:10, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1, # male or symptom off
    TRUE ~ .
  )))

################################ IMPUTATION ####################################
### Prepare the df and merge with og df
abcd_qcd <- abcd_qcd %>% dplyr::select(-c('pubertal_stage','matdep','income'))
template <- expand.grid(id = unique(abcd_qcd$id), time = 1:8)
prepped_df <- merge(template, abcd_qcd, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# reorder cols
prepped_df <- prepped_df[, c("id", "time", "sex", names(abcd_qcd)[3:8])]
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(`worthless`:`worry`, ~replace(., is.na(.), NA)))
# fill sex (have checked that individuals don't have conflicting sex)
prepped_df <- prepped_df %>%
  group_by(id) %>%
  fill(sex, .direction = "down") %>% fill(sex, .direction = "up")
# Check completed df
print(prepped_df)
### inspect missing data
plot_pattern(prepped_df, square = TRUE, rotate = TRUE)  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ABCD")
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

### make predictor matrix of symptoms only
predMat <- make.predictorMatrix(prepped_df)
predMat[,1] <- 0 # don't use ID
meth <- make.method(prepped_df)
meth[c("sex")] <- ""

#### impute
imputed <- mice(prepped_df, m=19, maxit=20, meth=meth, seed=123, predictorMatrix=predMat)
# extract and combine pooled datasets into long format
imputed_list <- lapply(1:19, function(i) complete(imputed, action = i))
all_imputed <- bind_rows(imputed_list)
# sanity check: count the number of duplicate rows based on 'id' and 'time' (11726*8 = 93808)
nrow(all_imputed %>%
       group_by(id, time, sex) %>%
       filter(n() > 1) %>%
       summarise(count = n()) %>%
       ungroup())

# group by 'id', 'time', 'sex' and the symptom columns, then count occurrences
combination_counts <- all_imputed %>%
  group_by(id,time,sex,worthless,anxious,guilty,self_conscious,unhappy,worry) %>%
  summarise(count = n(), .groups = 'drop')

# for each grouping of 'id', 'time', 'sex', select the combination with the highest count
abcd_mode_imputed <- combination_counts %>%
  group_by(id, time, sex) %>%
  slice_max(count, n = 1, with_ties = FALSE) %>%
  ungroup() %>% dplyr::select(-count)

# final dataset that pools the mode symptom combination from multiple imputation
head(abcd_mode_imputed) ## nrow = 93808

################################################################################

### sex stratification
boys <- abcd_mode_imputed %>% filter(sex==-1)
girls <- abcd_mode_imputed %>% filter(sex==1)

## [0,1] encoding
imputed_abcd_recode <- abcd_mode_imputed %>%
  mutate(across(4:9, ~case_when(
    . == 1 ~ 1,
    . == -1 ~ 0,
    TRUE ~ .
  )))

########################## MULTIGROUP ISING MODEL ##############################
# vars to use
vars <- names(abcd_mode_imputed)[4:ncol(abcd_mode_imputed)]

# Form saturated model and run [all params free]
model1 <- Ising(imputed_abcd_recode, vars = vars, groups = "time") %>% runmodel
# Prune-stepup to find a sparse model:
model1b <- model1 %>% prune(alpha = 0.05) %>%  stepup(alpha = 0.05)
# Equal networks (omega = network structure, edges btw nodes equal across time)
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
network_bpm <- getmatrix(best.model, "omega")[[1]]
graph_bpm <- qgraph(network_bpm, layout = 'spring', labels = vars, theme = 'colorblind')

# overlay centrality plot
centralityPlot(list(Wave1 = network_bpm),
               theme_bw=FALSE, scale = "z-scores", 
               include = c("Strength","Closeness","Betweenness"), 
               labels = vars, orderBy = 'Strength')

# extract temperature and calculate 95%CIs from standard errors of the beta parameter
temp_bpm <-  as.numeric(lapply(getmatrix(best.model, "beta"), 'mean'))
betas_est <- best.model@parameters$est[best.model@parameters$matrix == "beta"]
betas_se <- best.model@parameters$se[best.model@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_bpm + (z*betas_se)
lowerCI <- temp_bpm - (z*betas_se)

# extract external fields (information)
#fields_bpm <- lapply(getmatrix(best.model, 'tau'), 'mean')

# set labels
rownames(network_bpm) <- labels
colnames(network_bpm) <- labels

# plot heatmaps (2 ways)
heatmap(network_bpm, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA, main = paste("BPM items"))

heatmap.2(network_bpm, 
        symm = TRUE,
        col = viridis::plasma(100),Rowv = NA,
        trace = "none", density.info = "none")

## plotting
ages <- c(10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14)
#stages <- c('pre','early','mid','late','post')
#income <- c('<5k','5-12k','12-16k','16-25k','25-35k','35-50k','50-75k','75-100k','100-200k', '>200k')

# gather temp data for plotting
temp_data <- data.frame(Age = ages,
                        #Stage = stages,
                        #Income = income,
                        Temperature = 1/temp_bpm,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)

temp_data_girls$sex <- "Females"
temp_data_boys$sex <- "Males"
temp_data <- rbind(temp_data_boys, temp_data_girls)

# convert stage to factor levels to plot x axis in order
#temp_data$Stage <- factor(temp_data$Stage, levels = stages)
#temp_data$Income <- factor(temp_data$Income, levels = income)

# plot temp change with CIs
ggplot(temp_data, aes(x = Age, y = Temperature)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.55, 1.0) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),  
        plot.title = element_text(hjust = 0.5, size = 12), legend.position='none')  +
  scale_x_continuous(breaks = ages) +
  #scale_discrete_manual(breaks = income) + 
  scale_color_manual(values = c("Females" = "orange", "Males" = "blue"))

####################### BOOTSTRAPPING #######################

# read in results from eddie/datastore
setwd('/Volumes/igmm/GenScotDepression/users/poppy/abcd')
boot.abcd <- readRDS('bootstrap_results_abcd.rds')
print(boot.abcd)

std_errors <- c(NA,apply(boot.abcd$t, 2, sd))
z = qnorm(0.975)
upper_boot_CI <- temp_bpm + (z*std_errors)
lower_boot_CI <- temp_bpm - (z*std_errors)

ages <- c(10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14)
boot_ci_data <- data.frame(Age = ages,
                           Temperature = 1/temp_bpm,
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
  ylim(0.55, 1.0) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), 
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position='NULL') + 
  scale_x_continuous(breaks = ages) +
  scale_color_manual(values = c('Analytical 95% CIs' = 'purple', 'Bootstrapped 95% CIs' = 'red'), name=NULL) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "solid"))))



################## SENSITIVITY/SUPPLEMENTAL ANALYSES ########################

## histograms of overall depression score 
abcd_qcd$total <- rowSums(abcd_qcd[,3:8])

breaks <- seq(min(abcd_qcd$total, na.rm = TRUE), max(abcd_qcd$total, na.rm = TRUE), length.out = 8)

hist(abcd_qcd$total[abcd_qcd$time==2], main = 'Age 11', xlab = 'Overall depression', breaks=breaks)
mtext('(c)', 3, at = -46, padj = -4)

par(mar = c(6, 2, 6, 4))

hist(abcd_qcd$total[abcd_qcd$time==8], main = 'Age 14', xlab = 'Overall depression', breaks=breaks)

dev.off()

## distribution of sumscores
# Add a new col 'sumscore' 
imputed_bpm_recode <- imputed_bpm_recode %>%
  dplyr::mutate(sumscore = rowSums(dplyr::select(., worthless:worry)))

# plot histograms of sumscore across time (age)
ggplot(imputed_bpm_recode, aes(x = sumscore)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  facet_wrap(~time, ncol = 4) +
  labs(title = "Distribution of Sumscore for [0,1] Encoding by Time in ABCD Study",
       x = "Sumscore",
       y = "Frequency") + geom_vline(xintercept = 0, linetype = "dashed", color = "black")


#### other stratification groups

## puberty stratification
pre <- imputed_bpm %>% filter(pubertal_stage==1)
early <- imputed_bpm %>% filter(pubertal_stage==2)
mid <- imputed_bpm %>% filter(pubertal_stage==3)
late <- imputed_bpm %>% filter(pubertal_stage==4)
post <- imputed_bpm %>% filter(pubertal_stage==5)

## income stratification
income <- c('<5k','5-12k','12-16k','16-25k','25-35k','35-50k','50-75k','75-100k','100-200k', '>200k')
bracket1 <- imputed_bpm %>% filter(income==1)
bracket3 <- imputed_bpm %>% filter(income==3)
bracket5 <- imputed_bpm %>% filter(income==5)
bracket7 <- imputed_bpm %>% filter(income==7)
bracket9 <- imputed_bpm %>% filter(income==9)
