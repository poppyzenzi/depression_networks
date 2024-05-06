##############################
## ALSPAC network sensitivity
##############################
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
##############################

# Function to simulate Ising network data (2 datasets/time points)
simulate_data <- function(n, p, beta, dataset_id) {
  # Generate random binary data
  data <- matrix(sample(c(-1, 1), n * p, replace = TRUE), nrow = n)
  # Compute pairwise correlations based on Ising model
  cor_matrix <- diag(p)
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      cor_matrix[i, j] <- cor_matrix[j, i] <- mean(data[, i] * data[, j])
    }
  }
  # Apply inverse temperature parameter (beta)
  cor_matrix <- tanh(beta * cor_matrix)
  # Add dataset_id column
  data <- cbind(data, dataset_id)
  
  return(list(data = data, cor_matrix = cor_matrix))
}

# Gradually change characteristics of the simulated datasets
temp_values <- seq(1, 0.5, by = -0.02) # from 1 to 0.5, increments of 0.02
beta_values <- 1/temp_values

combined_data <- NULL
first_dataset <- simulate_data(n = 100, p = 13, beta = beta_values[1], dataset_id = 1)
combined_data <- rbind(combined_data, cbind(first_dataset$data, beta = beta_values[1]))

for (i in seq_along(beta_values[-1])) {
  beta <- beta_values[i + 1]
  dataset <- simulate_data(n = 100, p = 13, beta = beta, dataset_id = 2)
  dataset$cor_matrix <- tanh(beta * first_dataset$cor_matrix)  # Apply beta to the correlation matrix
  combined_data <- rbind(combined_data, cbind(dataset$data, beta = beta))
}

combined_data <- as.data.frame(combined_data)
colnames(combined_data) <- c(paste0("var", 1:13), "dataset_id", "beta")
head(combined_data)

#####################################################################

# Variables to use:
vars <- names(combined_data)[1:13]

# Form saturated model and run [all params free]
model1 <- Ising(combined_data, vars = vars, groups = "beta")
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

#extract temperature
temp_test <- as.numeric(lapply(getmatrix(model2, "beta"), 'mean'))
# calculate 95%CIs from standard errors of the beta parameter
betas_est <- model2@parameters$est[model2b@parameters$matrix == "beta"]
betas_se <- model2@parameters$se[model2b@parameters$matrix == "beta"]
z = qnorm(0.975)
upperCI <- temp_test + (z*betas_se)
lowerCI <- temp_test - (z*betas_se)

## plot temperature change
temp_data <- data.frame(beta = beta_values,
                        Temperature = 1/temp_test,
                        UpperCI = 1/upperCI,
                        LowerCI = 1/lowerCI)

ggplot(temp_data, aes(x = beta, y = Temperature)) +
  geom_point(color = "blue", shape = 16) +
  geom_line() +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) +
  labs(x = "Age", y = "Temperature", title = "Temperature change")




#################################################
################### symptom skew ################

setwd('/Volumes/igmm/GenScotDepression/users/poppy/alspac')
smfq_dat <- read.table('smfq_symptoms_wide.txt', check.names = FALSE)
labels <- c("unhappy", "anhedonia", "apathetic", "restless", "worthless",
            "tearful", "distracted", "self_loathing", "guilty", "isolated", 
            "unloved", "inadequate", "incompetent")
colnames(smfq_dat) <- c('id', 'time', unlist(labels))
## filter for only valid symptom scores btw 0-3
smfq_symptoms <- smfq_dat %>%
  filter(if_all(.cols = all_of(names(smfq_dat)[3:ncol(smfq_dat)]), ~ . >= 0 & . <= 3))

data_subset <- smfq_symptoms[, 3:15]
par(mfrow=c(4, 4)) 

for (i in 1:ncol(data_subset)) {
  hist(data_subset[, i], breaks = 0:3, main = colnames(data_subset)[i], xlab = "Values", cex.main=0.9)
}

main_title <- "Histograms of SMFQ Scores in ALSPAC"
mtext(main_title, outer = TRUE, cex = 0.8, font = 2, line = -1)

par(mfrow=c(1, 1))

#################################################

