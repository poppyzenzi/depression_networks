#### ALSPAC SDQ temperature 

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
# ages: 7, 9.6, 12, 13, 16.5
## read in the data, add predictor variables and labels
setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
sdq_qcd <- read.table('alspac_sdq_wide_qcd.txt', check.names = FALSE)

labels <- c("malaise", "worries", "unhappy", "anxiety","fears", "solitary",
            "friends","liked","bullied", "adultoriented")

emo.labels <- c("malaise", "worries", "unhappy", "anxiety","fears")

colnames(sdq_qcd)[4:13] <- unlist(labels)

sdq.emo <- sdq_qcd[1:8]

########## IMPUTATION ##############

### Prepare the df and merge with og df
template <- expand.grid(id = unique(sdq.emo$id), time = 2:5) # exclude t1 (age7)
prepped_df <- merge(template, sdq.emo, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(malaise:fears, ~replace(., is.na(.), NA)))
# fill sex (have checked that individuals don't have conflicitng sex)
prepped_df <- prepped_df %>%
  group_by(id) %>%
  fill(sex, .direction = "down") %>% fill(sex, .direction = "up")
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
predMat[,c(1:3)] <- 0
meth <- make.method(prepped_df)
meth[c("sex")] <- ""

#### impute
imputed_sdq <- mice(prepped_df,m=5,maxit=40,meth=meth, predictorMatrix=predMat)
summary(imputed_sdq)
imputed_df_sdq <- complete(imputed_sdq,1)

write.table(imputed_df_sdq, 'alspac_imputed_ising_sdq.txt', col.names=TRUE)

# Variables to use:
vars <- names(imputed_df_sdq)[4:ncol(imputed_df_sdq)]
boys <- imputed_df_sdq %>% filter(sex==1)
girls <- imputed_df_sdq %>% filter(sex==2)

# Form saturated model and run [all params free]
model1 <- Ising(imputed_df_sdq, vars = vars, groups = "time")
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
best.model <- model3

### CHANGE TO SDQ
#extract and plot network
all_network_sdq <- getmatrix(best.model, "omega")
network_sdq <- getmatrix(best.model, "omega")[[1]]
graph_sdq <- qgraph(network_sdq, layout = 'spring', labels = vars, theme = 'colorblind')

## overlay centrality plot
centralityPlot(list(Wave1 = all_network_sdq[[1]]),
theme_bw=FALSE, scale = "z-scores", 
include = c("Strength","Closeness","Betweenness"), 
labels = vars, orderBy = 'Strength')

#extract temperature
temp_sdq <- as.numeric(lapply(getmatrix(best.model, "beta"), 'mean'))
# calculate 95%CIs from standard errors of the beta parameter
betas_est <- best.model@parameters$est[best.model@parameters$matrix == "beta"]
betas_se <- best.model@parameters$se[best.model@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_sdq + (z*betas_se)
lowerCI <- temp_sdq - (z*betas_se)
## plot temperature change
ages <- c(10,12,13,17)
temp_data <- data.frame(Age = ages,
                        Temperature = 1/temp_sdq,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)


ggplot(temp_data, aes(x = Age, y = Temperature)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.8, 1.1) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), 
        plot.title = element_text(hjust = 0.5, size = 12)) + 
  scale_x_continuous(breaks = ages)

rownames(network_sdq) <- labels[1:5]
colnames(network_sdq) <- labels[1:5]

heatmap(network_sdq, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA,
        main = paste("SDQ items"))
