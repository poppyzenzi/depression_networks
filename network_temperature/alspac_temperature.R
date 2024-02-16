#########################
## ALSPAC network temp
#########################
library (foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library (qgraph)

#########################
#Data

setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
#smfq_qcd <- read.table('smfq_symptoms_qcd.txt', check.names = FALSE)

# all variables
smfq_qcd <- read.table('network_all_vars.txt', check.names = FALSE)

labels <- c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
            "tearful", "distracted", "self-loathing", "guilty", "isolated", 
            "unloved", "inadequate", "incompetent", "sex",
            "maternal depression","bullying", "child trauma", "sleep", "income", "MDDPRS")


colnames(smfq_qcd) <- c('id', 'time', unlist(labels))

#recode variables so that each variable is binary with +1 and -1 
smfq_qcd <- smfq_qcd %>%
  mutate(across(3:19, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1,
    TRUE ~ .
  )))

# remove second wave
filtered <- smfq_qcd[smfq_qcd$time != 2, ]
smfq_qcd <- filtered

### fit psychonetrics model

#fit an Ising model with increasing constraints representing their hypotheses to this longitudinal assessment 
# We investigate the impact on the fit of the model of 
# 1) constraining edges between nodes to be equal across time points
# 2) constraining the external fields to be equal across time points 
# 3) constraining the temperature (the entropy of the system) to be equal across time points. 
# Additionally, we tested whether a dense network (all nodes are connected) or 
#  a sparse network (at least some edges are absent) fits the data best

# Variables to use:
vars <- names(smfq_qcd)[3:19]

# Form saturated model and run [all params free]
model1 <- Ising(smfq_qcd, vars = vars, groups = "time")
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
network_smfq <- getmatrix(model2, "omega")[[1]]
graph_smfq <- qgraph(network_smfq, layout = 'spring', labels = vars, theme = 'colorblind')

#extract temperature 
temp_smfq <-  as.numeric(lapply(getmatrix(model2, "beta"), 'mean'))

#extract external fields
fields_smfq <- lapply(getmatrix(model2, 'tau'), 'mean')


### plotting

pdf('alspac_network_temp.pdf', 7, 5)
layout(matrix(c(1,1,2,2,2,2,2,3,3,3,3,3,
                4,4,4,4,5,5,5,5,6,6,6,6), 2, 12, byrow = TRUE))

par(mar = rep(0,4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
                   "tearful", "distracted", "self-loathing", "guilty", "isolated", 
                   "unloved", "inadequate", "incompetent"), 
       title = 'SMFQ items', col = 'darkorange', pch = 19,
       cex = 0.8, bty = 'n')

# a label at 3rd margin (top) at 0.1 along and right justification (2)
mtext('(a)', 3, at = .01, padj = 2)

qgraph(network_smfq, layout = 'spring', 
       groups = list(1:13), palette = 'colorblind',
       legend = FALSE, theme = 'colorblind',
       labels = c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
                  "tearful", "distracted", "self-loathing", "guilty", "isolated", 
                  "unloved", "inadequate", "incompetent"), vsize = 12)

par(mar = rep(2,4), cex.main = 0.8)

plot(1/temp_smfq, bty = 'n', xlab = 'Age', ylab = 'Temperature', xaxt = 'n', yaxt = 'n', 
     #ylim = c(.85, 1), 
     type = 'b', main = 'Change in network temperature')
axis(1, c(seq(1, 3, 1)), c('11', '14', '17'))
axis(2, c(seq(.8, 2, .05)))
mtext('(b)', 3, at = .32, padj = -2)

par(mar = c(6, 4, 6, 2))

# histograms of overall depression score 
smfq_qcd$total <- rowSums(smfq_qcd[,3:15])

hist(smfq_qcd$total[smfq_qcd$time==1], main = 'Age 11', xlab = 'Overall depression')
mtext('(c)', 3, at = -46, padj = -4)

par(mar = c(6, 3, 6, 3))

#hist(smfq_qcd$total[smfq_qcd$time==2], main = 'Age 13', xlab = 'Overall depression',  yaxt = 'n', ylab = '')
hist(smfq_qcd$total[smfq_qcd$time==3], main = 'Age 14', xlab = 'Overall depression',  yaxt = 'n', ylab = '')

par(mar = c(6, 2, 6, 4))

hist(smfq_qcd$total[smfq_qcd$time==4], main = 'Age 17', xlab = 'Overall depression',  yaxt = 'n', ylab = '')

dev.off()