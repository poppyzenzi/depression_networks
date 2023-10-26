#################################################################
# "Cross-Lagged Network Models"
# using total bpm score at each time point as the indicator 
#################################################################

library(glmnet)
library(qgraph)
library(lavaan)
library(dplyr)

# load and clean data 
setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/") 
# load wide symptom data, all waves, 3 point or binary
all_waves <- read.table(file = "bpm_sum_wide.txt", sep = "", header = TRUE) %>% .[,c(-1)]
all_waves[all_waves == -9999] <- NA 

#################################################################

#text for legends (to keep track of what each variable is)

#9 I feel worthless or inferior
#11 I am too fearful or anxious	
#12 I feel too guilty	
#13 I am self-conscious or easily embarrassed	
#18 I am unhappy, sad, or depressed	
#19 I worry a lot

longlegend <- c("9: I feel worthless or inferior", 
                "11: I am too fearful or anxious", 
                "12: I feel too guilty", 
                "13: I am self-conscious or easily embarrassed", 
                "18: I am unhappy, sad, or depressed", 
                "19: I worry a lot")

brieflegend <- c("9: worthless", 
                 "11: fearful", 
                 "12: guilty", 
                 "13: self-conscious", 
                 "18: unhappy", 
                 "19: worry")
#################################################################
#################################################################

# remove cases with missing data (glmnet doesn't allow them)
df <- all_waves[complete.cases(all_waves),]
nrow(df) #3588 cases remain

#################################################################

# estimate CLPM with latent variables ###########################

# configural = same structure across waves but allows loadings, intercepts and residual variances to vary
# metric (weak inv) = constrains factor loadings to equality across waves
# scalar = intercepts constrained to equality across waves (e.g. £100 inflation)

configural <- '
        # define latent factors and factor loadings
          dep.1 =~ NA*t1 + t2 + t3 + t4 + t5 + t6 + t7 + t8

        # Variances and Covariance
          dep.1 ~~ dep.1

        # residual variances - allow them to be estimated freely
          t1 ~~ t1
          t2 ~~ t2
          t3 ~~ t3
          t4 ~~ t4
          t5 ~~ t5
          t6 ~~ t6
          t7 ~~ t7
          t8 ~~ t8
'

fit_configural <- cfa(configural, data = df, mimic = "mplus")
summary(fit_configural, fit.measures = TRUE)
semPaths(fit_configural, what="est", sizeLat = 7, sizeMan = 5, edge.label.cex = .8)


metric <- '# define latent factors and factor loadings
          dep.1 =~ NA*bpm9_y1 + lambda1*bpm11_y1 + lambda2*bpm12_y1 + lambda3*bpm13_y1 + lambda4*bpm18_y1 + lambda5*bpm19_y1
          dep.2 =~ NA*bpm9_y2 + lambda1*bpm11_y2 + lambda2*bpm12_y2 + lambda3*bpm13_y2 + lambda4*bpm18_y2 + lambda5*bpm19_y2
          
        # Intercepts
          bpm9_y1 ~ i1*1
          bpm11_y1 ~ 1
          bpm12_y1 ~ 1
          bpm13_y1 ~ 1
          bpm18_y1 ~ 1
          bpm19_y1 ~ 1
          bpm9_y2 ~ i1*1
          bpm11_y2 ~ 1
          bpm12_y2 ~ 1
          bpm13_y2 ~ 1
          bpm18_y2 ~ 1
          bpm19_y2 ~ 1

          # Unique Variances
          bpm9_y1 ~~ bpm9_y1
          bpm11_y1 ~~ bpm11_y1
          bpm12_y1 ~~ bpm12_y1
          bpm13_y1 ~~ bpm13_y1
          bpm18_y1 ~~ bpm19_y1
          bpm9_y2 ~~ bpm9_y2
          bpm11_y2 ~~ bpm11_y2
          bpm12_y2 ~~ bpm12_y2
          bpm13_y2 ~~ bpm13_y2
          bpm18_y2 ~~ bpm19_y2
          
          # latent variable means
          dep.1 ~ 0*1
          dep.2 ~ 1
          
          # Latent Variable Variances and Covariance
          dep.1 ~~ 1*dep.1
          dep.2 ~~ dep.2
          dep.1 ~~ dep.2'

fit_metric <- cfa(metric, data = df, mimic = "mplus")
summary(fit_metric, fit.measures = TRUE)

semPaths(fit_metric, what="est",
         sizeLat = 7, sizeMan = 5, edge.label.cex = .8)


### CLPM model ####
config.model <- paste(CLPM.mod, sep=' ', collapse=NULL)
CLPM.fit <- sem(CLPM.mod, data = df, estimator = "MLM", std.lv = TRUE)

CLPM.fit <- sem(CLPM.mod, data = df, estimator = "MLM", std.lv = TRUE)
summary(CLPM.fit, standardized = TRUE,fit = TRUE)
#################################################################