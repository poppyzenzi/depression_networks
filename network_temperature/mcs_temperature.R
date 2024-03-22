#########################
## MCS network temp
#########################
library (foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library (qgraph)
library(mice)
library(ggmice)
library(VIM)

#########################
#Data

setwd('/Volumes/igmm/GenScotDepression/users/poppy/mcs/symptom_data')
mcs_qcd <- read.table('mcs_sdq_sym_long.txt', check.names = FALSE)

labels <- c("malaise", "worries", "unhappy", "anxiety","fears", "solitary",
            "friends*","liked*","bullied", "adult-oriented")

colnames(mcs_qcd) <- c('id', 'time', unlist(labels))

#recode variables so that each variable is binary with +1 and -1 
mcs_qcd <- mcs_qcd %>%
  mutate(across(3:12, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1,
    TRUE ~ .
  )))


########## IMPUTATION ##############

### Prepare the df and merge with og df
template <- expand.grid(id = unique(mcs_qcd$id), time = 4:7)
prepped_df <- merge(template, mcs_qcd, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(`malaise`:`adult-oriented`, ~replace(., is.na(.), NA)))
# reorder cols
prepped_df <- prepped_df[, c("id", "time", names(mcs_qcd)[3:ncol(mcs_qcd)])]
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
predMat[,c(1:2)] <- 0
meth <- make.method(prepped_df)

#### impute
imputed <- mice(prepped_df,m=3,maxit=50,meth=meth, predictorMatrix=predMat)
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
vars <- names(mcs_qcd)[3:12]

# Form saturated model and run [all params free]
model1 <- Ising(mcs_qcd, vars = vars, groups = "time", estimator = 'ML')
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

#extract and plot network
# "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
network_sdq <- getmatrix(model2, "omega")[[1]]
graph_sdq <- qgraph(network_sdq, layout = 'spring', labels = vars, theme = 'colorblind')

#extract temperature 
temp_sdq <-  as.numeric(lapply(getmatrix(model2, "beta"), 'mean'))

#extract external fields
fields_sdq <- lapply(getmatrix(model2, 'tau'), 'mean')


rownames(network_sdq) <- labels
colnames(network_sdq) <- labels

heatmap(network_sdq, 
        symm = TRUE,
        col = viridis::plasma(100),
        Rowv = NA)

### plotting

pdf('mcs_network_temp.pdf', 7, 5)
layout(matrix(c(1,1,2,2,2,2,2,3,3,3,3,3,
                4,4,4,4,5,5,5,5,6,6,6,6), 2, 12, byrow = TRUE))

par(mar = rep(0,4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', labels, 
       title = 'SDQ items', col = 'darkorange', pch = 19,
       cex = 0.8, bty = 'n')

# a label at 3rd margin (top) at 0.1 along and right justification (2)
mtext('(a)', 3, at = .01, padj = 2)

qgraph(network_sdq, layout = 'spring', 
       groups = list(1:10), palette = 'classic',
       legend = FALSE, theme = 'colorblind',
       labels = labels, vsize = 12)

par(mar = rep(2,4), cex.main = 0.8)

plot(1/temp_sdq, bty = 'n', xlab = 'Age', ylab = 'Temperature', xaxt = 'n', yaxt = 'n', 
     ylim = c(.85, 1), 
     type = 'b', main = 'Change in network temperature')
axis(1, c(seq(1, 4, 1)), c('7', '11', '14', '17'))
axis(2, c(seq(.8, 2, .05)))
mtext('(b)', 3, at = .32, padj = -2)

par(mar = c(6, 4, 6, 2))

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