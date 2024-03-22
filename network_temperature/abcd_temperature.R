#########################
## ABCD network temp
#########################
library (foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library (qgraph)

#########################
#Data

setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")
abcd_qcd <- read.table('abcd_bpm_sym_long.txt', check.names = FALSE)
labels <- c("worthless","anxious","guilty","self_conscious","unhappy","worry")

colnames(abcd_qcd) <- c('id', 'time', unlist(labels))

#recode variables so that each variable is binary with +1 and -1 
abcd_qcd <- abcd_qcd %>%
  mutate(across(3:8, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1,
    TRUE ~ .
  )))

# remove second wave
#filtered <- abcd_qcd[abcd_qcd$time %in% c(2, 8), ]
#abcd_qcd <- filtered

########## IMPUTATION ##############

### Prepare the df and merge with og df
template <- expand.grid(id = unique(abcd_qcd$id), time = 1:8)
prepped_df <- merge(template, abcd_qcd, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(`worthless`:`worry`, ~replace(., is.na(.), NA)))
# reorder cols
prepped_df <- prepped_df[, c("id", "time", names(abcd_qcd)[3:ncol(abcd_qcd)])]
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
predMat[,c(1:2)] <- 0
meth <- make.method(prepped_df)

#### impute
imputed <- mice(prepped_df,m=5,maxit=50, meth=meth,seed=500, predictorMatrix=predMat)
summary(imputed)
imputed_df <- complete(imputed,1)

### fit psychonetrics model

#fit an Ising model with increasing constraints representing their hypotheses to this longitudinal assessment 
# We investigate the impact on the fit of the model of 
# 1) constraining edges between nodes to be equal across time points
# 2) constraining the external fields to be equal across time points 
# 3) constraining the temperature (the entropy of the system) to be equal across time points. 
# Additionally, we tested whether a dense network (all nodes are connected) or 
#  a sparse network (at least some edges are absent) fits the data best

# Variables to use:
vars <- names(abcd_qcd)[3:8]

# Form saturated model and run [all params free]
model1 <- Ising(abcd_qcd, vars = vars, groups = "time")
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
network_bpm <- getmatrix(model2, "omega")[[1]]
graph_bpm <- qgraph(network_bpm, layout = 'spring', labels = vars, theme = 'colorblind')

#extract temperature 
temp_bpm <-  as.numeric(lapply(getmatrix(model2, "beta"), 'mean'))

#extract external fields
fields_bpm <- lapply(getmatrix(model2, 'tau'), 'mean')


rownames(network_bpm) <- labels
colnames(network_bpm) <- labels

heatmap(network_bpm, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA)

heatmap.2(network_bpm, 
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

qgraph(network_bpm, layout = 'spring', 
       groups = list(1:6), palette = 'colorblind',
       legend = FALSE, theme = 'colorblind',
       labels = labels, vsize = 12)

par(mar = rep(2,4), cex.main = 0.8)

plot(1/temp_bpm, bty = 'n', xlab = 'Age', ylab = 'Temperature', xaxt = 'n', yaxt = 'n', 
     ylim = c(.7, 1), 
     type = 'b', main = 'Change in network temperature')
axis(1, c(seq(1, 8, 1)), c('10.5', '11','11.5','12','12.5','13','13.5','14'))
axis(2, c(seq(.7, 2, .05)))
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