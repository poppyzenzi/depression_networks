#########################
## ABCD network temp [CBCL version]
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
#Data


setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")
abcd_qcd <- read.table('cbcl_symptoms_binary.txt', check.names = FALSE)
labels <- c('anhedonia','tearful','selfharm','eating','worthlessness',
            'guilt','fatigue','oversleeps','undersleeps','suicidalideation',
            'sleep','lethargic','unhappy')

# 0 - not true, 1 - true
#There is very little they enjoy
#Cries a lot
#Deliberately harms self or attempts suicide
#Doesnt eat well
#Feels worthless or inferior
#Feels too guilty
#Overtired without good reason
#Sleeps less than most kids
#Sleeps more than most kids during day and/or night
#Talks about killing self
#Trouble sleeping
#Underactive, slow moving, or lacks energy
#Unhappy, sad, or depressed

colnames(abcd_qcd) <- c('id', 'time', unlist(labels))

## demographics
demog <- read.table('/Volumes/igmm/GenScotDepression/users/poppy/abcd/abcd_envi_vars.txt')
sexdf <- demog %>% filter(eventname==0) %>% select(src_subject_id, demo_sex_v2)
names(sexdf) <- c('id','sex')
abcd_qcd <- merge(abcd_qcd, sexdf, by='id', all.x=TRUE)

#recode variables so that each variable is binary with +1 and -1 
abcd_qcd <- abcd_qcd %>%
  mutate(across(3:16, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1, # male
    TRUE ~ .
  )))

########## IMPUTATION ##############
### Prepare the df and merge with og df
template <- expand.grid(id = unique(abcd_qcd$id), time = c(0,2,4,6,8)) # give everyone 5 observations
prepped_df <- merge(template, abcd_qcd, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# reorder cols
prepped_df <- prepped_df[, c("id", "time", "sex", names(abcd_qcd)[3:15])]
# Fill missing sex values by id
prepped_df <- prepped_df %>%
  arrange(id, time) %>%
  group_by(id) %>%
  fill(sex, .direction = "downup") %>%
  ungroup()
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(`anhedonia`:`unhappy`, ~replace(., is.na(.), NA)))
# Check completed df
print(prepped_df)

### inspect 
plot_pattern(prepped_df, square = TRUE, rotate = TRUE)  + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ABCD")
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

### make predictor matrix
# we only want to impute the symptoms vars for now
predMat <- make.predictorMatrix(prepped_df)
predMat[,c(1:3)] <- 0
meth <- make.method(prepped_df)

#### impute
cbcl_imputed <- mice(prepped_df,m=5,maxit=50, meth=meth,seed=500, predictorMatrix=predMat)
summary(cbcl_imputed)
cbcl_imputed_df <- complete(cbcl_imputed,1)

### sex stratification
boys <- cbcl_imputed_df %>% filter(sex==-1)
girls <- cbcl_imputed_df %>% filter(sex==1)

### fit psychonetrics model

#fit an Ising model with increasing constraints representing their hypotheses to this longitudinal assessment 
# We investigate the impact on the fit of the model of 
# 1) constraining edges between nodes to be equal across time points
# 2) constraining the external fields to be equal across time points 
# 3) constraining the temperature (the entropy of the system) to be equal across time points. 
# Additionally, we tested whether a dense network (all nodes are connected) or 
#  a sparse network (at least some edges are absent) fits the data best

# Variables to use:
vars <- names(cbcl_imputed_df)[4:16]

# Form saturated model and run [all params free]
model1 <- Ising(girls, vars = vars, groups = "time")
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

#extract and plot network
# "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
network_cbcl <- getmatrix(model2, "omega")[[1]]
graph_cbcl <- qgraph(network_cbcl, layout = 'spring', labels = vars, theme = 'colorblind')

#extract temperature 
temp_cbcl <-  as.numeric(lapply(getmatrix(model2, "beta"), 'mean'))
# calculate 95%CIs from standard errors of the beta parameter
betas_est <- model2@parameters$est[model2b@parameters$matrix == "beta"]
betas_se <- model2@parameters$se[model2b@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_cbcl + (z*betas_se)
lowerCI <- temp_cbcl - (z*betas_se)

#extract external fields (information)
fields_cbcl <- lapply(getmatrix(model2, 'tau'), 'mean')

rownames(network_cbcl) <- labels
colnames(network_cbcl) <- labels

heatmap(network_cbcl, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA, main = paste("CBCL items"))

heatmap.2(network_cbcl, 
          symm = TRUE,
          col = viridis::plasma(100),Rowv = NA,
          trace = "none", density.info = "none")

### plotting

pdf('mcs_network_temp.pdf', 7, 5)
layout(matrix(c(1,2,2,2,2,3,3,3,3,
                4,4,4,4,4,5,5,5,5,5), 2, 9, byrow = TRUE))

par(mar = rep(0,4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', labels, 
       title = 'BPM items', col = 'darkorange', pch = 19,
       cex = 0.8, bty = 'n')

# a label at 3rd margin (top) at 0.1 along and right justification (2)
mtext('(a)', 3, at = .01, padj = 2)

qgraph(network_cbcl, layout = 'spring', 
       groups = list(1:6), palette = 'colorblind',
       legend = FALSE, theme = 'colorblind',
       labels = labels, vsize = 12)

par(mar = rep(2,4), cex.main = 0.8)

## plot temperature change

ages <- c(10, 11, 12, 13, 14)

# plot estimates
temp_data <- data.frame(Age = ages,
                        Temperature = 1/temp_cbcl,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)

# Plot estimates and confidence intervals using ggplot2
ggplot(temp_data, aes(x = Age, y = Temperature)) +
  geom_point(color = "blue", shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  ylim(0.85, 1.1) +
  labs(x = "Age", y = "Temperature", title = "Temperature change") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 11),  # Adjust size of x-axis tick labels
        axis.text.y = element_text(size = 11),  # Adjust size of y-axis tick labels
        plot.title = element_text(hjust = 0.5, size = 12)) +  # Center and adjust size of title
  scale_x_continuous(breaks = ages)

mtext('(b)', 3, at = .32, padj = -2)
par(mar = c(6, 4, 6, 2))

# histograms of overall depression score 
abcd_qcd$total <- rowSums(abcd_qcd[,3:8])

breaks <- seq(min(abcd_qcd$total, na.rm = TRUE), max(abcd_qcd$total, na.rm = TRUE), length.out = 8)

hist(abcd_qcd$total[abcd_qcd$time==2], main = 'Age 11', xlab = 'Overall depression', breaks=breaks)
mtext('(c)', 3, at = -46, padj = -4)

par(mar = c(6, 2, 6, 4))

hist(abcd_qcd$total[abcd_qcd$time==8], main = 'Age 14', xlab = 'Overall depression', breaks=breaks)

dev.off()