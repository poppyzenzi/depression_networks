#########################
## ALSPAC network temp
#########################
library (foreign)
library(bootnet)
library(psychonetrics)
library(dplyr)
library (qgraph)
library(mice)
library(ggmice)

#########################
#Data
setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
#smfq_qcd <- read.table('smfq_symptoms_qcd.txt', check.names = FALSE)

# all variables
smfq_qcd <- read.table('network_all_vars.txt', check.names = FALSE) %>%
  select(1:16,22:24)

labels <- c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
            "tearful", "distracted", "self_loathing", "guilty", "isolated", 
            "unloved", "inadequate", "incompetent", "sex",
             "mood", "psychotic", "neurodev")

colnames(smfq_qcd) <- c('id', 'time', unlist(labels))

#recode variables so that each variable is binary with +1 and -1 
smfq_qcd <- smfq_qcd %>%
  mutate(across(3:15, ~case_when(
    . == 1 ~ 1,
    . == 0 ~ -1,
    TRUE ~ .
  )))

########## IMPUTATION ##############

# remove non symptom vars
symptoms_only <- smfq_qcd %>% select(1:15)
### Prepare the df and merge with og df
template <- expand.grid(id = unique(symptoms_only$id), time = 1:7)
prepped_df <- merge(template, symptoms_only, by = c("id", "time"), all.x = TRUE)
# sort df by id and time
prepped_df <- prepped_df %>% arrange(id, time)
# NA for non-time invariant columns
prepped_df <- prepped_df %>%
  mutate(across(unhappy:incompetent, ~replace(., is.na(.), NA)))
# fill PRS for time invariant columns
#prepped_df <- prepped_df %>%
#  group_by(id) %>%
#  fill(mood:neurodev)
# reorder cols
prepped_df <- prepped_df[, c("id", "time", names(symptoms_only)[3:ncol(symptoms_only)])]
# Check completed df
print(prepped_df)

### inspect 
plot_pattern(prepped_df, square = TRUE, rotate = TRUE)
aggr_plot <- aggr(prepped_df, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(prepped_df), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

### make predictor matrix
# we only want to impute the symptoms vars for now
predMat <- make.predictorMatrix(prepped_df)
predMat[,c(1:2)] <- 0
meth <- make.method(prepped_df)

#### impute
imputed <- mice(prepped_df,m=5,maxit=50,meth=meth,seed=500, predictorMatrix=predMat)
summary(imputed)
imputed_df <- complete(imputed,1)

########## options to stratify the data ###########

# stratify by sex
girls <- subset(smfq_qcd, sex==1)
boys <- subset(smfq_qcd, sex==-1)

# genetically straitfy by quintiles
quintiles <- quantile(smfq_qcd$MDDPRS, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1))
labels <- c("very low", "low", "middle", "high", "very high")
smfq_qcd$quintile <- cut(smfq_qcd$MDDPRS, breaks = quintiles, labels = labels, include.lowest = TRUE)

# get distribution and frequency
hist(smfq_qcd$mood)
count <- table(smfq_qcd$quintile)
count

# get mean and SD of PRS in each group 
risk_group_types <- unique(smfq_qcd$quintile)  # Get quintiles

very_high <- subset(smfq_qcd, quintile=='very high')
very_low <- subset(smfq_qcd, quintile=='very low')
low <- subset(smfq_qcd, quintile=='low')
high <- subset(smfq_qcd, quintile=='high')
med <- subset(smfq_qcd, quintile=='middle')
high_quart <- subset(smfq_qcd, HighQuartile==1) # nunique = 1552

# dataframe for summary statistics
summary_statistics <- data.frame(Category = character(),
                                 Mean = numeric(),
                                 SD = numeric(),
                                 Count = integer(),  # Change to integer type
                                 stringsAsFactors = FALSE)
# Loop through each risk group type
for (risk_group_type in risk_group_types) {
  mean_value <- mean(smfq_qcd$MDDPRS[smfq_qcd$risk_group == risk_group_type], na.rm = TRUE)
  sd_value <- sd(smfq_qcd$MDDPRS[smfq_qcd$risk_group == risk_group_type], na.rm = TRUE)
  count_value <- count[risk_group_type]
  summary_statistics <- rbind(summary_statistics, 
                              data.frame(Category = risk_group_type, 
                                         Mean = mean_value, 
                                         SD = sd_value, 
                                         Count = count_value))
}

summary_statistics


## highest quartile method
quartiles <- quantile(smfq_qcd$MDDPRS, probs = c(0.8))  # 75th percentile (highest quartile)
smfq_qcd$HighQuartile <- ifelse(smfq_qcd$MDDPRS >= quartiles, 1, 0)
head(smfq_qcd)

######### subset complete cases only ######### 

# group df by 'id' and count time points per id
id_counts <- prepped_df %>%
  group_by(id) %>%
  summarise(n_time_points = n_distinct(time))
# filter to include indiviudals with all 7 tps
complete_ids <- id_counts %>%
  filter(n_time_points == 7) %>%
  pull(id)
# subset df
complete_df <- prepped_df %>%
  filter(id %in% complete_ids)
# print subsetted 
head(complete_df)
length(unique(complete_df$id)) # 1009

######### fit psychonetrics Ising model ##########

#fit an Ising model with increasing constraints representing their hypotheses to this longitudinal assessment 
# We investigate the impact on the fit of the model of 
# 1) constraining edges between nodes to be equal across time points
# 2) constraining the external fields to be equal across time points 
# 3) constraining the temperature (the entropy of the system) to be equal across time points. 
# Additionally, we tested whether a dense network (all nodes are connected) or 
#  a sparse network (at least some edges are absent) fits the data best

# Variables to use:
vars <- names(prepped_df)[3:15]

# Form saturated model and run [all params free]
model1 <- Ising(prepped_df, vars = vars, groups = "time", estimator = 'ML')

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

#extract and plot network
network_smfq <- getmatrix(model2, "omega")[[1]]
graph_smfq <- qgraph(network_smfq, layout = 'spring', labels = vars, theme = 'colorblind')

#extract temperature 
temp_smfq <-  as.numeric(lapply(getmatrix(model2, "beta"), 'mean'))

#extract external fields
fields_smfq <- lapply(getmatrix(model2b, 'tau'), 'mean')


### plotting

pdf('alspac_network_temp.pdf', 7, 5)
layout(matrix(c(1,1,2,2,2,2,2,3,3,3,3,3,
                4,4,4,4,5,5,5,5,6,6,6,6), 2, 12, byrow = TRUE))

par(mar = rep(0,4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend('center', c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
                   "tearful", "distracted", "self-loathing", "guilty", "isolated", 
                   "unloved", "inadequate", "incompetent"), 
       title = 'SMFQ items', pch = 19,
       cex = 0.8, bty = 'n')

# a label at 3rd margin (top) at 0.1 along and right justification (2)
mtext('(a)', 3, at = .01, padj = 2)

# palette (for nodes): "rainbow"(default), "colorblind","pastel","gray","R"and"ggplot2"
# theme (for edges): "classic", "colorblind", "gray", "Hollywood", "Borkulo", "gimme", "TeamFortress", "Reddit", "Leuven" or "Fried".

qgraph(network_smfq, layout = 'spring', 
       groups = list(1:13), 
       posCol='slateblue', negCol='red', palette='colorblind',
       legend = FALSE,
       labels = c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
                  "tearful", "distracted", "self-loathing", "guilty", "isolated", 
                  "unloved", "inadequate", "incompetent"), vsize = 10)

par(mar = rep(2,4), cex.main = 0.8)

plot(1/temp_smfq, bty = 'n', xlab = 'Age', ylab = 'Temperature', xaxt = 'n', yaxt = 'n', 
     ylim = c(.85, 1.05), 
     type = 'b', main = 'Change in network temperature')
axis(1, c(seq(1, 7, 1)), c('11', '13','14', '17','18','19','22'))
axis(2, c(seq(.7, 2, .05)))
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


### Compare networks of data sets using NCT ###
## Networks can be compared by either (1) feeding the data directly into NCT (whereby 
## you need to specify arguments such as "gamma" and "binary.data") or (2) by using 
## estimateNetwork() (bootnet package) and feeding that output into NCT. For the latter 
## option, we refer to the help file of  estimateNetwork() for its usage. Below, both 
## options are illustrated. We recommend using estimateNetwork(), since this function 
## has implemented many network estimation methods.

## gamma = 0 (in estimateNetwork this hyperparameter is called "tuning"; to illustrate 
# how to specify a different value than the default)
## iterations (it) set to 10 to save time
## Note: Low number of iterations can give unreliable results; should be 1000 at least

## Testing whether there are differences in the three aspects that are validated 
# (network invariance, global strength, edge weight)
## 2 edges are tested here: between variable 1 and 2, and between 3 and 6 (can be 
# "list(c(2,1),c(6,3))" as well)

## (1) Feeding data directly into NCT, plot results of global strength invariance test
set.seed(123)
NCT_a <- NCT(subset(smfq_qcd, time==1)[,3:15], subset(smfq_qcd, time==2)[,3:15], gamma=0, it=50, binary.data = TRUE, 
             test.edges=TRUE, edges=list(c(1,2),c(3,6)))
summary(NCT_a)
plot(NCT_a, what="strength")

NCT_b <- NCT(subset(smfq_qcd, time==2)[,3:15], subset(smfq_qcd, time==3)[,3:15], gamma=0, it=50, binary.data = TRUE, 
             test.edges=TRUE, edges=list(c(1,2),c(3,6)))
summary(NCT_b)
plot(NCT_b, what="strength")

NCT_c <- NCT(subset(smfq_qcd, time==3)[,3:15], subset(smfq_qcd, time==4)[,3:15], gamma=0, it=50, binary.data = TRUE, 
             test.edges=TRUE, edges=list(c(1,2),c(3,6)))
summary(NCT_c)
plot(NCT_c, what="strength")

###### code from Jonas psychonetrics example #####

# Extract beta:
beta <-  unlist(getmatrix(model2b, "beta"))

# Standard errors::
SEs <- model2b@parameters$se[model2b@parameters$matrix == "beta"]

# Make a data frame:
df <- data.frame(
  temperature = 1/beta,
  group = names(beta),
  lower = 1 / (beta-qnorm(0.975) * SEs),
  upper = 1 / (beta+qnorm(0.975) * SEs),
  stringsAsFactors = FALSE
)

# Some extra values:
df$fixed <- is.na(df$lower)
df$group <- factor(df$group)

# Create the plot:
library("ggplot2")
g <- ggplot(df,aes(x=as.numeric(group), y = temperature, ymin = lower, ymax = upper)) + 
  geom_line() + 
  geom_errorbar(width = 0.05) + 
  geom_point(cex = 5, colour = "black") +
  geom_point(aes(colour = fixed), cex = 4) +  theme_bw()  +
  xlab("") + ylab(expression(paste("Temperature (",1/beta,")"))) + 
  scale_x_continuous(breaks =  1:7, labels = levels(df$group), expand = c(0.1,0.1)) + 
  scale_y_continuous(expand = c(0,.1), limits = c(0,1)) + 
  theme( panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+ 
  ggtitle(expression(paste("Model 6: 2 groups; ",bold(Omega)," sparse & equal; ",bold(tau)," equal; ",beta," free"))) + 
  scale_colour_manual(values = c("black","white")) + 
  theme(legend.position = "none")

# Plot:
print(g)

# Make labels:
labels <- c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
            "tearful", "distracted", "self-loathing", "guilty", "isolated", 
            "unloved", "inadequate", "incompetent")

# Extract network structure and thresholds:
network <- getmatrix(model2b, "omega")[[1]]
thresholds <- getmatrix(model2b, "tau")[[1]]

# Scale thresholds for colors:
scaledthresh <- as.vector(thresholds / (2*max(abs(thresholds))))

# Make colors:
cols <- ifelse(scaledthresh < 0, "red", "darkblue")
cols[scaledthresh>0] <- qgraph:::Fade(cols[scaledthresh>0],alpha = scaledthresh[scaledthresh>0], "white")
cols[scaledthresh<0] <- qgraph:::Fade(cols[scaledthresh<0],alpha = abs(scaledthresh)[scaledthresh<0], "white")

# Plot network and save to file:
qgraph(network, layout = "spring", labels = labels,
       shape = "rectangle", vsize = 15, vsize2 = 8,
       theme = "colorblind", color = cols,
       cut = 0.5, repulsion = 0.9)

