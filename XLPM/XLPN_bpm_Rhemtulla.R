#################################################################
#Code to accompany Rhemtulla, van Bork, & Cramer (in preparation) 
# "Cross-Lagged Network Models"
#################################################################
library(glmnet)
library(qgraph)
library(lavaan)

# load and clean data ###########################################
setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/") 
# load wide symptom data, all waves, 3 point or binary
all_waves <- read.table(file = "bpm_symptoms_binary.txt", sep = "", header = FALSE)

colnames(all_waves) <- c(
  "bpm9_y1", "bpm9_y2", "bpm9_y3", "bpm9_y4", "bpm9_y5",
  "bpm9_y6", "bpm9_y7", "bpm9_y8", "bpm11_y1", "bpm11_y2",
  "bpm11_y3", "bpm11_y4", "bpm11_y5", "bpm11_y6", "bpm11_y7",
  "bpm11_y8", "bpm12_y1", "bpm12_y2", "bpm12_y3", "bpm12_y4",
  "bpm12_y5", "bpm12_y6", "bpm12_y7", "bpm12_y8", "bpm13_y1",
  "bpm13_y2", "bpm13_y3", "bpm13_y4", "bpm13_y5", "bpm13_y6",
  "bpm13_y7", "bpm13_y8", "bpm18_y1", "bpm18_y2", "bpm18_y3",
  "bpm18_y4", "bpm18_y5", "bpm18_y6", "bpm18_y7", "bpm18_y8",
  "bpm19_y1", "bpm19_y2", "bpm19_y3", "bpm19_y4", "bpm19_y5",
  "bpm19_y6", "bpm19_y7", "bpm19_y8"
)

#change missingness indicator from -9999 (in original data) to NA
all_waves[all_waves == -9999] <- NA 

#################################################################

# collapse highly similar variable if desired##########################

#reorder variables to put all bpm wave 1 variables together, etc.: 
all_waves <- all_waves[,c(1, 9, 17, 25, 33, 41, 2, 10, 18, 26, 34, 42, 3, 11,
                          19, 27, 35, 43, 4, 12, 20, 28, 36, 44, 5, 13, 21, 29,
                          37, 45, 6, 14, 22, 30, 38, 46, 7, 15, 23, 31, 39, 47, 
                          8, 16, 24, 32, 40, 48)]

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

# save data and legends (not necessary)###########################
save(all_waves, longlegend, brieflegend, file = "all_waves.Rda")
load(file = "all_waves.Rda") #to load the data and legends without running all of the above code
#################################################################

# select waves we want, here 2 and 4
df <- all_waves[,c(7:12,19:24)]

# remove cases with missing data (glmnet doesn't allow them)
df <- df[complete.cases(df),]
nrow(df) #9756 cases remain

colnames(df) <- c('bpm9_y1', 'bpm11_y1','bpm12_y1','bpm13_y1', 'bpm18_y1','bpm19_y1',
                  'bpm9_y2','bpm11_y2','bpm12_y2','bpm13_y2','bpm18_y2','bpm19_y2')

# estimate network paths using glmnet##########################################

k <- 6 #k is the number of variables at each time point, 6 symptoms
adjMat <- matrix(0, k, k) #set up empty matrix of coefficients
rsquarelist <- rep(0, k)
#the following code loops the regularized regression through each
#variable at time 2 as a DV, regressing on each variable at time 1
#as predictors. It assumes that the first k columns in your data
#correspond to the variables at time 1 and the next k correspond
#to the variables at time 2. These variables must be in the same 
#order at each time point.
#please see the paper for technical details

for (i in 1:k){
  set.seed(100)
  lassoreg <- cv.glmnet(as.matrix(df[,1:k]), df[,(k+i)], 
                        family = "gaussian", alpha = 1, standardize=TRUE)
  lambda <- lassoreg$lambda.min 
  #  rsquare[i] <- lassoreg$glmnet.fit$dev.ratio[which(lassoreg$lambda == lambda)]
  adjMat[1:k,i] <- coef(lassoreg, s = lambda, exact = FALSE)[2:(k+1)]
}
#################################################################

# plot the CLPN ##################################################

#create a 'group' vector that says which variables belong to which construct:
#groups <- c(rep("depression", 6)) 

#create a vector of variable labels for the figure
labels <- c("9","11","12","13","18","19")

nodecol = "royalblue"
#jpeg("CLPN_ARCLpaths.jpg", width=8, height=4, res=800, units="in")  #(uncomment this code to save plot as jpeg)
layout(mat = matrix(c(1, 2), 1, 2), widths = c(1.5, 1)) #to plot multiple graphs in one figure
qgraph(adjMat, labels = labels, legend = FALSE, colors = nodecol)
plot.new()
par(mar = rep(0, 4))
legend(x = "center", inset = c(0, 0), bty = "n", legend = brieflegend,  
       col = nodecol, y.intersp = 2, pch = 19) #add legend to identify variables
#dev.off() #(uncomment this code to save plot as jpeg)


# curved arrows on each variable are the AR paths
# we can remove those to see the CL paths more clearly: 

adjMat2 <- adjMat
diag(adjMat2) <- 0 #set the AR paths to 0 for plotting purposes: 

#jpeg("CLPN_CLpaths.jpg", width=8, height=4, res=800, units="in")  #(uncomment this code to save plot as jpeg)
layout(mat = matrix(c(1, 2), 1, 2), widths = c(1.5, 1)) #to plot multiple graphs in one figure
qgraph(adjMat2, labels = labels, legend = FALSE, colors = nodecol)
plot.new()
par(mar = rep(0, 4))
legend(x = "center", inset = c(0, 0), bty = "n", legend = brieflegend,  
       col = nodecol, y.intersp = 2, pch = 19) #add legend to identify variables
#dev.off() #(uncomment this code to save plot as jpeg)

#################################################################
# compute and plot predictability and influence using lavaan #################

# predictability of node A is extent we can predict A by its neighbours N(A)
# nodewise predictability = r-squared based on regressions of 
# each variable on (1) all others at previous timepoint (CL paths only)
# and (2) based on only other construct(s) at previous timepoint 

pred0 <- rep(0, k) #set up empty vectors to hold predictability coefficients
pred1 <- rep(0, k) 
#pred2 <- rep(0, k) # only one group of latent factor 

for (i in 1:k) {  # loop through each variable
  name.i <- names(df)[k + 1] #what is the name of the DV
  include.mod0 <- names(df)[1:k]  # Include self at prev time point
  include.mod1 <- names(df)[1:k][-i] # Exclude self at previous time point
  #include.mod2 <- names(df)[1:k]   # exclude as only one 'group' of latent factors
  betas.mod0 <- adjMat[1:k, i]  # Fix regression coefficients to the ones obtained by glmnet
  betas.mod1 <- adjMat[1:k, i][-i]
  #betas.mod2 <- adjMat[1:k, i] # excluded as only one 'group' of latent factors 
  # Construct regression formulas: DV ~ (betas*predictors)
  mod0 <- paste(name.i, " ~ ", paste(paste(betas.mod0,'*', include.mod0, sep = ""), collapse = "+")) 
  mod1 <- paste(name.i, " ~ ", paste(paste(betas.mod1,'*', include.mod1, sep = ""), collapse = "+")) 
  #mod2 <- paste(name.i, " ~ ", paste(paste(betas.mod2,'*', include.mod2, sep = ""), collapse = "+")) 
  # sample covariance matrix, covariance matrix, n observations, fit SEM
  fit0 <- sem(mod0, sample.cov = cor(df), sample.nobs = nrow(df)) #fit model with subset of predictors to get r-square
  pred0[i] <- inspect(fit0, "rsquare") # inspect gets specific info from lavaan object, get SEM model and rsq val
  fit1 <- sem(mod1, sample.cov = cor(df), sample.nobs = nrow(df)) #fit model with subset of predictors to get r-square
  pred1[i] <- inspect(fit1, "rsquare")
  #fit2 <- sem(mod2, sample.cov = cor(df), sample.nobs = nrow(df)) #fit model with subset of predictors to get r-square
  #pred2[i] <- inspect(fit2, "rsquare")
}

pred0
pred1
#pred2

#nodewise influence: sum squared outgoing regression paths/node
infl0 <- apply(adjMat^2, 1, sum) # all outgoing nodes including AR
infl1 <- apply(adjMat2^2, 1, sum) # all outgoing nodes excluding AR

#adjMat3 <- adjMat2 #adjMat3 sets paths within the same construct to 0
#adjMat3[groups[row(adjMat3)] == groups[col(adjMat3)]] <- 0
#infl2 <- apply(adjMat3^2, 1, sum) #rows = T1 predictors; columns = T2 DVs

infl0
infl1
#infl2



#create a plot of predictability and influence: 
#jpeg("PredictabilityInfluencePlot.jpg", width=8, height=6, res=800, units="in")  #(uncomment this code to save plot as jpeg)
layout(mat = matrix(c(1, 2, 3, 3), 2, 2), widths = c(1, 1)) #to plot multiple graphs in one figure
par(mar = c(2, 4, 0, 0), las = 1) #set margins
barplot(rev(pred1), horiz = TRUE, names.arg = rev(labels),
        col = viridis(k), ylab = "predictability")
barplot(rev(infl1), horiz = TRUE, names.arg = rev(labels), 
        col = viridis(k), ylab = "influence")
plot.new()
legend(x = "center", inset = c(0, 0), bty = "n", fill=rev(viridis(k)),
       legend = brieflegend, cex = 1.2) #add legend to identify variables
#dev.off() #(uncomment this code to save plot as jpeg)
#################################################################
# make corr plot of symptoms at two waves

library(corrplot)
corrplot(cor(df[,1:12]), tl.col='black', tl.cex = .75)

# estimate CLPM with latent variables ###########################

# configural = same structure across waves but allows loadings, intercepts and residual variances to vary
# metric (weak inv) = constrains factor loadings to equality across waves
# scalar = intercepts constrained to equality across waves (e.g. £100 inflation)

configural <- '
        # define latent factors and factor loadings
          dep.1 =~ NA*bpm9_y1 + bpm11_y1 + bpm12_y1 + bpm13_y1 + bpm18_y1 + bpm19_y1
          dep.2 =~ NA*bpm9_y2 + bpm11_y2 + bpm12_y2 + bpm13_y2 + bpm18_y2 + bpm19_y2
      
        # Variances and Covariance
          dep.1 ~~ dep.1
          dep.2 ~~ dep.2
          dep.1 ~~ dep.2
          
        # residual variances - allow them to be estimated freely
          bpm9_y1 ~~ bpm9_y1
          bpm11_y1 ~~ bpm11_y1
          bpm12_y1 ~~ bpm12_y1
          bpm13_y1 ~~ bpm13_y1
          bpm18_y1 ~~ bpm18_y1
          bpm19_y1 ~~ bpm19_y1
          bpm9_y2 ~~ bpm9_y2
          bpm11_y2 ~~ bpm11_y2
          bpm12_y2 ~~ bpm12_y2
          bpm13_y2 ~~ bpm13_y2
          bpm18_y2 ~~ bpm18_y2
          bpm19_y2 ~~ bpm19_y2

          bpm9_y1 ~~ bpm9_y2
          bpm11_y1 ~~ bpm11_y2
          bpm12_y1 ~~ bpm12_y2
          bpm13_y1 ~~ bpm13_y2
          bpm18_y1 ~~ bpm18_y2
          bpm19_y1 ~~ bpm19_y2
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