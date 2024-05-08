#########################
## ABCD network temp
#########################
library (foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library (qgraph)
library(VIM)
library(mice)
library(ggmice)

#########################
# data (symptoms at t1,2,3,4,5,6,7,8)
setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")
abcd_qcd <- read.table('abcd_bpm_sym_long.txt', check.names = FALSE)
labels <- c("worthless","anxious","guilty","self_conscious","unhappy","worry")

colnames(abcd_qcd) <- c('id', 'time', unlist(labels))

# demographics
demog <- read.table('/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_envi_vars.txt')
env.vars <- demog %>% filter(eventname==0) %>% select(src_subject_id, demo_sex_v2, fam_history_q6d_depression, demo_comb_income_v2)
names(env.vars) <- c('id','sex','matdep','income')
abcd_qcd <- merge(abcd_qcd, env.vars, by='id', all.x=TRUE)

# puberty (at t0,2,4,6,8, t0 is lost as no symptoms)
puberty <- read.table('../predictors/puberty.vars.txt') %>% select(src_subject_id, eventname, pds_y_ss)
names(puberty) <- c('id','time','pubertal_stage')
abcd_qcd <- merge(abcd_qcd, puberty, by=c('id','time'), all.x=TRUE)

#recode variables so that each variable is binary with +1 and -1 
abcd_qcd <- abcd_qcd %>%
  mutate(across(3:10, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1, # male or symptom off
    TRUE ~ .
  )))
 
########## IMPUTATION ##############
### Prepare the df and merge with og df
template <- expand.grid(id = unique(abcd_qcd$id), time = 1:8)
prepped_df <- merge(template, abcd_qcd, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# reorder cols
prepped_df <- prepped_df[, c("id", "time", "sex", "pubertal_stage", "matdep", "income", names(abcd_qcd)[3:8])]
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(`worthless`:`worry`, ~replace(., is.na(.), NA)))
# Check completed df
print(prepped_df)
### inspect 
plot_pattern(prepped_df, square = TRUE, rotate = TRUE)  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ABCD")
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

### make predictor matrix
# we only want to impute the symptoms vars for now, not pubertal stage etc.
predMat <- make.predictorMatrix(prepped_df)
predMat[,c(1:6)] <- 0
meth <- make.method(prepped_df)
meth[c("sex", "pubertal_stage","matdep","income")] <- ""

#### impute
imputed <- mice(prepped_df,m=5,maxit=50, meth=meth, predictorMatrix=predMat)
summary(imputed)
imputed_bpm <- complete(imputed,1)

### sex stratification
boys <- imputed_bpm %>% filter(sex==-1)
girls <- imputed_bpm %>% filter(sex==1)

### fit psychonetrics model

#fit an Ising model with increasing constraints representing their hypotheses to this longitudinal assessment 
# We investigate the impact on the fit of the model of 
# 1) constraining edges between nodes to be equal across time points
# 2) constraining the external fields to be equal across time points 
# 3) constraining the temperature (the entropy of the system) to be equal across time points. 
# Additionally, we tested whether a dense network (all nodes are connected) or 
#  a sparse network (at least some edges are absent) fits the data best

# Variables to use:
vars <- names(imputed_bpm)[7:ncol(imputed_bpm)]
# drop NA rows for puberty info
imputed_puberty <- na.omit(imputed_df, cols = "pubertal_stage")
boys <- imputed_puberty %>% filter(sex==-1)
girls <- imputed_puberty %>% filter(sex==1)

# Form saturated model and run [all params free]
model1 <- Ising(na.omit(imputed_bpm, cols='income'), vars = vars, groups = "income")
model1 <- model1 %>% runmodel
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

best.model <- model3

#extract and plot network
# "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
network_bpm <- getmatrix(best.model, "omega")[[1]]
graph_bpm <- qgraph(network_bpm, layout = 'spring', labels = vars, theme = 'colorblind')

## overlay centrality plot
centralityPlot(list(Wave1 = network_bpm),
               theme_bw=FALSE, scale = "z-scores", 
               include = c("Strength","Closeness","Betweenness"), 
               labels = vars, orderBy = 'Strength')

#extract temperature 
temp_bpm <-  as.numeric(lapply(getmatrix(best.model, "beta"), 'mean'))
# calculate 95%CIs from standard errors of the beta parameter
betas_est <- best.model@parameters$est[best.model@parameters$matrix == "beta"]
betas_se <- best.model@parameters$se[best.model@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_bpm + (z*betas_se)
lowerCI <- temp_bpm - (z*betas_se)

#extract external fields (information)
fields_bpm <- lapply(getmatrix(best.model, 'tau'), 'mean')

rownames(network_bpm) <- labels
colnames(network_bpm) <- labels

heatmap(network_bpm, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA, main = paste("BPM items"))

heatmap.2(network_bpm, 
        symm = TRUE,
        col = viridis::plasma(100),Rowv = NA,
        trace = "none", density.info = "none")

### plotting

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', labels, 
       title = 'BPM items', col = 'darkorange', pch = 19,
       cex = 0.8, bty = 'n')

qgraph(network_bpm, layout = 'spring', 
       groups = list(1:6), palette = 'colorblind',
       legend = FALSE, theme = 'colorblind',
       labels = labels, vsize = 12)

## plot temperature change
ages <- c(10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14)
stages <- c('pre','early','mid','late','post')
income <- c('<5k','5-12k','12-16k','16-25k','25-35k','35-50k','50-75k','75-100k','100-200k', '>200k')

# plot estimates
temp_data <- data.frame(#Age = ages,
                        #Stage = stages,
                        Income = income,
                        Temperature = 1/temp_bpm,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)

temp_data_girls$sex <- "Females"
temp_data_boys$sex <- "Males"
temp_data <- rbind(temp_data_boys, temp_data_girls)

# convert stage to factor levels to plot x axis in order
temp_data$Stage <- factor(temp_data$Stage, levels = stages)
temp_data$Income <- factor(temp_data$Income, levels = income)

# Plot estimates and confidence intervals using ggplot2
ggplot(temp_data, aes(x = Income, y = Temperature)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.8, 1.3) +
  labs(x = "Income bracket ($)", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),  
        plot.title = element_text(hjust = 0.5, size = 12))  +
  scale_discrete_manual(breaks = income) + 
  scale_color_manual(values = c("Females" = "orange", "Males" = "blue"))



ggplot(temp_data, aes(x = income, y = Temperature)) +
  geom_point(shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.75, 1.0) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11), 
        plot.title = element_text(hjust = 0.5, size = 12)) + 
  scale_x_continuous(breaks = ages) + 
  scale_color_manual(values = c("Mat dep" = "green", "No mat dep" = "deeppink"))

# histograms of overall depression score 
abcd_qcd$total <- rowSums(abcd_qcd[,3:8])

breaks <- seq(min(abcd_qcd$total, na.rm = TRUE), max(abcd_qcd$total, na.rm = TRUE), length.out = 8)

hist(abcd_qcd$total[abcd_qcd$time==2], main = 'Age 11', xlab = 'Overall depression', breaks=breaks)
mtext('(c)', 3, at = -46, padj = -4)

par(mar = c(6, 2, 6, 4))

hist(abcd_qcd$total[abcd_qcd$time==8], main = 'Age 14', xlab = 'Overall depression', breaks=breaks)

dev.off()