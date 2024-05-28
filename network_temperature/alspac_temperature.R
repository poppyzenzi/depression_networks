#############################################################################
############################# ALSPAC network temp ##############################
############################################################################# 
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

# recode variables [-1,1]
smfq_qcd <- smfq_qcd %>%
  mutate(across(3:15, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1,
    TRUE ~ .
  )))

################################ IMPUTATION ####################################
smfq_qcd <- smfq_qcd %>% dplyr::select(-c('maternal_depression','bullying','child_trauma','sleep','income'))
### Prepare the df and merge with og df
template <- expand.grid(id = unique(smfq_qcd$id), time = 1:6) # remove t7 at imputation
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
prepped_df <- prepped_df[, c("id", "time", "sex", unlist(labels))]
# Check completed df
print(prepped_df)

### inspect 
plot_pattern(prepped_df, square = TRUE, rotate = TRUE) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ALSPAC")
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

### make predictor matrix for symptom vars
predMat <- make.predictorMatrix(prepped_df)
predMat[,1] <- 0
meth <- make.method(prepped_df)
meth[c("sex")] <- ""

#### impute (takes a while)
imputed_smfq <- mice(prepped_df,m=41,maxit=20,meth=meth,seed=123,predictorMatrix=predMat)
saveRDS(imputed_smfq, 'imputed_smfq_41.rds')

# extract and combine pooled datasets into long format
imputed_list_smfq <- lapply(1:41, function(i) complete(imputed_smfq, action = i))
all_imputed_smfq <- bind_rows(imputed_list_smfq)
# sanity check: count the number of duplicate rows based on 'id' and 'time' (9217*7 = 64519)
nrow(all_imputed_smfq %>%
       group_by(id, time, sex) %>%
       filter(n() > 1) %>%
       summarise(count = n()) %>%
       ungroup())

# group by 'id', 'time', 'sex' and the symptom columns, then count occurrences
combination_counts <- all_imputed_smfq %>%
  group_by(id,time,sex,unhappy,anhedonia,apathetic,restless,worthless,tearful,
           distracted,self_loathing,guilty,isolated,unloved,inadequate,incompetent) %>%
  summarise(count = n(), .groups = 'drop')

# for each grouping of 'id', 'time', 'sex', select the combination with the highest count
alspac_mode_imputed <- combination_counts %>%
  group_by(id, time, sex) %>%
  slice_max(count, n = 1, with_ties = FALSE) %>%
  ungroup() %>% dplyr::select(-count)

# final dataset that pools the mode symptom combination from multiple imputation
head(alspac_mode_imputed) ## nrow = 64519

write.table(alspac_mode_imputed, 'alspac_imputed_ising_mode.txt', col.names=TRUE)

########## options to stratify the data ###########

# stratify by sex
girls <- subset(alspac_mode_imputed, sex==1)
boys <- subset(alspac_mode_imputed, sex==0)

# sensitivity [0,1] encoding
# recode variables so that each symptoms is binary with +1 and -1 
imputed_smfq_recode <- alspac_mode_imputed %>%
  mutate(across(4:16, ~case_when(
    . == 1 ~ 1,
    . == -1 ~ 0,
    TRUE ~ .
  )))

########################## MULTIGROUP ISING MODEL ##############################

# vars
vars <- names(alspac_mode_imputed)[4:16]

# remove age 22 (t7)
no.t7.og <- smfq_qcd[smfq_qcd$time != 7, ] %>% as.data.frame() # non-imputed (complete cases)
no.t7.imp <- alspac_mode_imputed[alspac_mode_imputed$time != 7, ] %>% as.data.frame()
boys.no.t7 <- no.t7.imp %>% filter(sex==0)
girls.no.t7 <- no.t7.imp %>% filter(sex==1)

# Form saturated model and run [all params free] (estimator = 'ML')
model1 <- Ising(imputed_smfq_recode, vars = vars, groups = "time") %>% runmodel
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

#extract and plot network
all_network_smfq <- getmatrix(best.model, "omega")
network_smfq <- getmatrix(best.model, "omega")[[1]]
graph_smfq <- qgraph(network_smfq, layout = 'spring', labels = vars, theme = 'colorblind')

# centrality plot
centralityPlot(list(Wave1 = all_network_smfq[[1]]),
               theme_bw=FALSE, scale = "z-scores", 
               include = c("Strength","Closeness","Betweenness"), 
               labels = vars, orderBy = 'Strength')


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

#extract temperature and calculate 95%CIs from standard errors of the beta parameter
temp_smfq <- as.numeric(lapply(getmatrix(best.model, "beta"), 'mean'))
betas_est <- best.model@parameters$est[best.model@parameters$matrix == "beta"]
betas_se <- best.model@parameters$se[best.model@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_smfq + (z*betas_se)
lowerCI <- temp_smfq - (z*betas_se)

### plot temperature change
ages <- c(11,13,14,17,18,19)
temp_data <- data.frame(Age = ages,
                        Temperature = 1/temp_smfq,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)

temp_data_girls$sex <- "Females"
temp_data_boys$sex <- "Males"
temp_data <- rbind(temp_data_boys, temp_data_girls)

ggplot(temp_data, aes(x = Age, y = Temperature)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.6, 1.0) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), 
        plot.title = element_text(hjust = 0.5, size = 12), legend.position='none') + 
  scale_x_continuous(breaks = ages) + 
  scale_color_manual(values = c("Females" = "orange", "Males" = "blue"))


#extract external fields
fields_smfq <- lapply(getmatrix(model2, 'tau'), 'mean')


####################### BOOTSTRAPPING #######################

# read in results from eddie/datastore
setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
boot.alspac <- readRDS('bootstrap_results_alspac.rds')
print(boot.alspac)

std_errors <- c(NA,apply(boot.alspac$t, 2, sd))
z = qnorm(0.975)
upper_boot_CI <- temp_smfq + (z*std_errors)
lower_boot_CI <- temp_smfq - (z*std_errors)

ages <- c(11,13,14,17,18,19)
boot_ci_data <- data.frame(Age = ages,
                           Temperature = 1/temp_smfq,
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
  ylim(0.6, 1.0) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), 
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position='NULL') + 
  scale_x_continuous(breaks = ages) +
  scale_color_manual(values = c('Analytical 95% CIs' = 'purple', 'Bootstrapped 95% CIs' = 'red'), name=NULL) +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "solid"))))




################################################
################### plotting ################### 



# histograms of overall depression score 
imputed_df_smfq_smfq$total <- rowSums(imputed_df_smfq_smfq[,3:15])

hist(imputed_df_smfq_smfq$total[imputed_df_smfq_smfq$time==1], main = 'Age 11', xlab = 'Overall depression')
mtext('(c)', 3, at = -46, padj = -4)

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
