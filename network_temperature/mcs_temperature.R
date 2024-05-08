#########################
## MCS network temp
#########################
library(foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library(qgraph)
library(mice)
library(ggmice)
library(VIM)

#########################
#Data
setwd('/Volumes/igmm/GenScotDepression/users/poppy/mcs/symptom_data')
mcs_qcd <- read.table('mcs_sdq_sym_long.txt', check.names = FALSE)

labels <- c("malaise", "worries", "unhappy", "anxiety","fears", "solitary",
            "friends","liked","bullied", "adultoriented")

emo.labels <- c("malaise", "worries", "unhappy", "anxiety","fears")

colnames(mcs_qcd) <- c('id', 'sex', 'time', unlist(labels))

#recode variables so that each variable is binary with +1 and -1 
#friends and liked are reverse coded
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

# emotional scale only
mcs.emo <- mcs_qcd[1:8]

########## IMPUTATION ##############

### Prepare the df and merge with og df
template <- expand.grid(id = unique(mcs.emo$id), time = 5:7) # try 4:7 and 5:7?
prepped_df <- merge(template, mcs.emo, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(`malaise`:`fears`, ~replace(., is.na(.), NA)))
# fill sex (have checked that individuals don't have conflicitng sex)
prepped_df <- prepped_df %>%
  group_by(id) %>%
  fill(sex, .direction = "down") %>% fill(sex, .direction = "up")
## a lot of missing sex data?
# Check completed df
print(prepped_df)

### inspect 
plot_pattern(prepped_df) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("MCS")

aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

### make predictor matrix
# we only want to impute the symptoms vars for now
predMat <- make.predictorMatrix(prepped_df)
predMat[,c(1:3)] <- 0
meth <- make.method(prepped_df)
meth[c("sex")] <- ""

#### impute
# if this doesn't work, usually an issue with variable labels, change to plain text and should run 
imputed <- mice(prepped_df,m=4,maxit=50,meth=meth, predictorMatrix=predMat)
summary(imputed)
imputed_df_mcs <- complete(imputed,1)

### sex stratification
boys <- imputed_df_mcs %>% filter(sex==0)
girls <- imputed_df_mcs %>% filter(sex==1)

### fit psychonetrics model

#fit an Ising model with increasing constraints representing their hypotheses to this longitudinal assessment 
# We investigate the impact on the fit of the model of 
# 1) constraining edges between nodes to be equal across time points
# 2) constraining the external fields to be equal across time points 
# 3) constraining the temperature (the entropy of the system) to be equal across time points. 
# Additionally, we tested whether a dense network (all nodes are connected) or 
#  a sparse network (at least some edges are absent) fits the data best

# Variables to use:
vars <- names(imputed_df_mcs[4:ncol(imputed_df_mcs)])

# remove age 7 (t4 from t4-7)
no.t4 <- imputed_df_mcs[imputed_df_mcs$time != 4, ]
boys.no.t4 <- no.t4 %>% filter(sex==0)
girls.no.t4 <- no.t4 %>% filter(sex==1)

# Form saturated model and run [all params free]
model1 <- Ising(girls, vars = vars, groups = "time")
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
# "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
network_sdq <- getmatrix(best.model, "omega")[[1]]
graph_sdq <- qgraph(network_sdq, layout = 'spring', labels = vars, theme = 'colorblind')

## overlay centrality plot
all_network_sdq <- getmatrix(best.model, "omega")
centralityPlot(list(Wave1 = all_network_sdq[[1]]),
theme_bw=FALSE, scale = "z-scores", 
include = c("Strength","Closeness","Betweenness"), 
labels = vars, orderBy = 'Strength')

# extract temperature 
temp_sdq <-  as.numeric(lapply(getmatrix(best.model, "beta"), 'mean'))
# calculate 95%CIs from standard errors of the beta parameter
betas_est <- best.model@parameters$est[best.model@parameters$matrix == "beta"]
betas_se <- best.model@parameters$se[best.model@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_sdq + (z*betas_se)
lowerCI <- temp_sdq - (z*betas_se)


# extract external fields
fields_sdq <- lapply(getmatrix(best.model, 'tau'), 'mean')

# heatmaps
rownames(network_sdq) <- emo.labels
colnames(network_sdq) <- emo.labels

heatmap(network_sdq, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA,
        main = paste("SDQ items"))


### plotting

qgraph(network_sdq, layout = 'spring', 
       legend = FALSE, theme = 'colorblind',
       labels = labels, vsize = 12)

## plot temperature change
ages <- c(11,14,17)
# plot estimates
temp_data_girls <- data.frame(Age = ages,
                        Temperature = 1/temp_sdq,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)

temp_data_girls$sex <- "Females"
temp_data_boys$sex <- "Males"
temp_data <- rbind(temp_data_boys, temp_data_girls)

# Plot estimates and confidence intervals using ggplot2
ggplot(temp_data, aes(x = Age, y = Temperature, color=sex)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.75, 1.1) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),  # Adjust size of x-axis tick labels
        axis.text.y = element_text(size = 11),  # Adjust size of y-axis tick labels
        plot.title = element_text(hjust = 0.5, size = 12)) +  # Center and adjust size of title
  scale_x_continuous(breaks = ages) +
  scale_color_manual(values = c("Females" = "orange", "Males" = "blue"))



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