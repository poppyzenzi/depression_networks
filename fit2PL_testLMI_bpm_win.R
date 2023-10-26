' testing mirt package for LMI with dichotomous data'

library(mirt)
library(glmnet)
library(qgraph)
library(lavaan)
library(sirt)

setwd("/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/") 
#setwd("C:/Users/s2421111/OneDrive - University of Edinburgh/Edinburgh/networks/mplus_ri_clpm")

all_waves <- read.table(file = "bpm_symptoms_binary.txt", sep = "", header = FALSE)

colnames(all_waves) <- c(
  "bpm09_y1", "bpm09_y2", "bpm09_y3", "bpm09_y4", "bpm09_y5",
  "bpm09_y6", "bpm09_y7", "bpm09_y8", "bpm11_y1", "bpm11_y2",
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

# collapse highly similar variable if desired##########################

#reorder variables to put all bpm wave 1 variables together, etc.: 
all_waves <- all_waves[,c(1, 9, 17, 25, 33, 41, 2, 10, 18, 26, 34, 42, 3, 11,
                          19, 27, 35, 43, 4, 12, 20, 28, 36, 44, 5, 13, 21, 29,
                          37, 45, 6, 14, 22, 30, 38, 46, 7, 15, 23, 31, 39, 47, 
                          8, 16, 24, 32, 40, 48)]


# select waves we want, here 2, 4, 6, 8 
df <- all_waves[,c(7:12,
                   19:24
                   #31:36
                   )]

# remove cases with missing data (glmnet doesn't allow them)
df <- df[complete.cases(df),]
nrow(df) #9756 cases remain


dat <- df
I <- ncol(dat)

###

# define Q-matrix
Q <- matrix(0,I,2)
Q[1:6,1] <- 1
Q[7:12,2] <- 1
#Q[13:18,3] <- 1
#Q[19:24,4] <- 1
rownames(Q) <- colnames(dat)
colnames(Q) <- c("T1","T2")

# vector with same items
itemnr <- as.numeric(substring( colnames(dat),4,5) )
# fix mean at T2 to zero
mu.fixed <- cbind( 2,0 )

# model 1: 2-dimensional 2PL model

# set variance at T2 to 1
variance.fixed <- cbind(2,2,1)

# M2a: rasch.mml2 (in sirt)
mod2a <- sirt::rasch.mml2(dat, Q=Q, est.b=itemnr, est.a=itemnr, mu.fixed=mu.fixed,
                          variance.fixed=variance.fixed, mmliter=100)
summary(mod2a)

#item_params_dim1 <- coef(mod2a, "est")$dim1
#item_params_dim2 <- coef(mod2a, "est")$dim2
#plot(item_params_dim1, item_params_dim2, xlab = "dim1", ylab="dim2")


#*************************************************
# Model 3: Concurrent calibration by assuming invariant item parameters
#*************************************************

 # use mirt for concurrent calibration
I <- ncol(dat)

# create user defined function for between item dimensionality 4PL model
name <- "4PLbw"
par <- c("low"=0,"upp"=1,"a"=1,"d"=0,"dimItem"=1)
est <- c(TRUE, TRUE,TRUE,TRUE,FALSE)
# item response function
irf <- function(par,Theta,ncat){
  low <- par[1]
  upp <- par[2]
  a <- par[3]
  d <- par[4]
  dimItem <- par[5]
  P1 <- low + ( upp - low ) * plogis( a*Theta[,dimItem] + d )
  cbind(1-P1, P1)
}


# create item response function
fourPLbetw <- mirt::createItem(name, par=par, est=est, P=irf)
head(dat)

# create mirt model (use variable names in mirt.model)
mirtsyn <- "
     T1=bpm09_y2,bpm11_y2,bpm12_y2,bpm13_y2,bpm18_y2,bpm19_y2
     T2=bpm09_y4,bpm11_y4,bpm12_y4,bpm13_y4,bpm18_y4,bpm19_y4
     COV=T1*T2,,T2*T2
     MEAN=T1
     CONSTRAIN=(bpm09_y2,bpm09_y4,d),(bpm11_y2,bpm11_y4,d),(bpm12_y2,bpm12_y4,d),
                (bpm13_y2,bpm13_y4,d),(bpm18_y2,bpm18_y4,d),(bpm19_y2,bpm19_y4,d),
                (bpm09_y2,bpm09_y4,a),(bpm11_y2,bpm11_y4,a),(bpm12_y2,bpm12_y4,a),
                (bpm13_y2,bpm13_y4,a),(bpm18_y2,bpm18_y4,a),(bpm19_y2,bpm19_y4,a),
        "

# create mirt model
mirtmodel <- mirt::mirt.model( mirtsyn, itemnames=colnames(dat) )
# define parameters to be estimated
mod3.pars <- mirt::mirt(dat, mirtmodel, rep( "4PLbw",I),
                        customItems=list("4PLbw"=fourPLbetw), pars="values")
# select dimensions
ind <- intersect( grep("T2",mod3.pars$item), which( mod3.pars$name=="dimItem" ) )
mod3.pars[ind,"value"] <- 2
# set item parameters low and upp to non-estimated
ind <- which( mod3.pars$name %in% c("low","upp") )
mod3.pars[ind,"est"] <- FALSE

# estimate 2PL model
mod3 <- mirt::mirt(dat, mirtmodel, itemtype=rep( "4PLbw",I),
                   customItems=list("4PLbw"=fourPLbetw), pars=mod3.pars, verbose=TRUE,
                   technical=list(NCYCLES=50)  )
mirt.wrapper.coef(mod3)


#####################################
# in lavaan
lavmodel <- "
             #**** T1
             F1=~ a1*bpm09_y2+a2*bpm11_y2+a3*bpm12_y2+a4*bpm13_y2+a5*bpm18_y2+a6*bpm19_y2
             bpm09_y2 | b1*t1 ; bpm11_y2 | b2*t1 ; bpm12_y2 | b3*t1 ; bpm13_y2 | b4*t1
             bpm18_y2 | b5*t1 ; bpm19_y2 | b6*t1
             F1 ~~ 1*F1
             #**** T2
             F2=~ a1*bpm09_y4+a2*bpm11_y4+a3*bpm12_y4+a4*bpm13_y4+a5*bpm18_y4+a6*bpm19_y4
             bpm09_y4 | b1*t1 ; bpm11_y4 | b2*t1 ; bpm12_y4 | b3*t1 ; bpm13_y4 | b4*t1
             bpm18_y4 | b5*t1 ; bpm19_y4 | b6*t1
             F2 ~~ NA*F2
             F2 ~ 1
             #*** covariance
             F1 ~~ F2
                "

mod3lav <- lavaan::cfa( data=dat, model=lavmodel,
                        std.lv=TRUE, ordered=colnames(dat), parameterization="theta")
summary(mod3lav, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)



#*************************************************
# Model 4: Linking with items of different item slope groups
#*************************************************

# dataset for T1
dat1 <- dat[, grep( "y2", colnames(dat) ) ]
colnames(dat1) <- gsub("y2","", colnames(dat1) )
# dataset for T2
dat2 <- dat[, grep( "y4", colnames(dat) ) ]
colnames(dat2) <- gsub("y4","", colnames(dat2) )

# 2PL model with slope groups T1
mod1 <- sirt::rasch.mml2( dat1, est.a=c( rep(1,2), rep(2,4) ) )
summary(mod1)

# 2PL model with slope groups T2
mod2 <- sirt::rasch.mml2( dat2, est.a=c( rep(1,4), rep(2,2) ) )
summary(mod2)

#------- Link 1: Haberman Linking
# collect item parameters
dfr1 <- data.frame( "study1", mod1$item$item, mod1$item$a, mod1$item$b )
dfr2 <- data.frame( "study2", mod2$item$item, mod2$item$a, mod2$item$b )
colnames(dfr2) <- colnames(dfr1) <- c("study", "item", "a", "b" )
itempars <- rbind( dfr1, dfr2 )
# Linking
link1 <- sirt::linking.haberman(itempars=itempars)
# R2 value close to 1 indicates strong invariance btw time points
# sqrtU2 values are root mean sq differences and should be small for invariance
# based on result: high R2 and low sqrtU2 and small DIF = high invariance, none excluded

#------- Link 2: Invariance alignment method
# create objects for invariance.alignment
nu <- rbind(c(mod1$item$thresh), c(mod2$item$thresh))
lambda <- rbind(c(mod1$item$a), c(mod2$item$a))
colnames(lambda) <- colnames(nu) <- c('bpm09','bpm11','bpm12','bpm13','bpm18','bpm19')
rownames(lambda) <- rownames(nu) <- c("T1", "T2")
# Linking
link2a <- sirt::invariance.alignment( lambda, nu )
summary(link2a)

