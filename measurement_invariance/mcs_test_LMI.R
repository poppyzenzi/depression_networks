#################################################################
# Millenium Cohort Study 
# testing for longitudinal measurement invariance
#################################################################

library(mirt)
library(glmnet)
library(qgraph)
library(lavaan)
library(sirt)

setwd("/Volumes/igmm/GenScotDepression/users/poppy") 

# read in binary symptom dat in wide format
all_waves <- read.table(file = "mcs/symptom_data/mcs_sdq_sym_wide.txt")

# drop ID col and rename symptom cols
all_waves <- all_waves[,c(-1)]
colnames(all_waves) <- c(
  "T5_1", "T5_2", "T5_3", "T5_4", "T5_5",
  "T5_6", "T5_7", "T5_8", "T5_9", "T5_10",
  "T6_1", "T6_2", "T6_3", "T6_4", "T6_5",
  "T6_6", "T6_7", "T6_8", "T6_9", "T6_10",
  "T7_1", "T7_2", "T7_3", "T7_4", "T7_5",
  "T7_6", "T7_7", "T7_8", "T7_9", "T7_10")

#change missingness indicator from -9999 (in original data) to NA
all_waves[all_waves == -9999] <- NA
all_waves[all_waves == -9] <- NA 


# select waves we want T5 T6 T7
df <- all_waves[,c(1:10,
                   21:30
                   #31:36
)]

# remove cases with missing data (glmnet doesn't allow them)
df <- df[complete.cases(df),]
nrow(df) #9191 cases remain

dat <- df
I <- ncol(dat)

###############

# define Q-matrix
Q <- matrix(0,I,2)
Q[1:10,1] <- 1
Q[11:20,2] <- 1
rownames(Q) <- colnames(dat)
colnames(Q) <- c("T1","T2")

# vector with same items, symptom/item numbers
itemnr <- as.numeric(substring(colnames(dat),4,5))
# fix mean at T2 to zero
mu.fixed <- cbind( 2,0 )

#*************************************************
#model 1: 2-dimensional 2PL model ##
#*************************************************

# set variance at T2 to 1
variance.fixed <- cbind(2,2,1)

# M2a: rasch.mml2 (in sirt)
mod2a <- sirt::rasch.mml2(dat, Q=Q, est.b=itemnr, est.a=itemnr, mu.fixed=mu.fixed,
                          variance.fixed=variance.fixed, mmliter=100)
summary(mod2a)

#item_params_dim1 <- coef(mod2a, "est")$dim1
#item_params_dim2 <- coef(mod2a, "est")$dim2
#plot(item_params_dim1, item_params_dim2, xlab = "dim1", ylab="dim2")

'
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
     T1=T5_1, T5_2, T5_3, T5_4, T5_5, T5_6, T5_7, T5_8, T5_9, T5_10
     T2=T6_1, T6_2, T6_3, T6_4, T6_5, T6_6, T6_7, T6_8, T6_9, T6_10
     COV=T1*T2,,T2*T2
     MEAN=T1
     CONSTRAIN=(T5_1,T6_1,d),(T5_2,T6_2,d),(T5_3,T6_3,d),(T5_4,T6_4,d),(T5_5,T6_5,d),
                (T5_6,T6_6,d),(T5_7,T6_7,d),(T5_8,T6_8,d),(T5_9,T6_9,d),(T5_10,T6_10,d),
                (T5_1,T6_1,a),(T5_2,T6_2,a),(T5_3,T6_3,a),(T5_4,T6_4,a),(T5_5,T6_5,a),
                (T5_6,T6_6,a),(T5_7,T6_7,a),(T5_8,T6_8,a),(T5_9,T6_9,a),(T5_10,T6_10,a),
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

# estimate 2PL model (slow)
mod3 <- mirt::mirt(dat, mirtmodel, itemtype=rep( "4PLbw",I),
                   customItems=list("4PLbw"=fourPLbetw), pars=mod3.pars, verbose=TRUE,
                   technical=list(NCYCLES=50)  )
mirt.wrapper.coef(mod3)


#####################################
# in lavaan
lavmodel <- "
             #**** T1
             F1=~ a1*T5_1+a2*T5_2+a3*T5_3+a4*T5_4+a5*T5_5+a6*T5_6+a7*T5_7+a8*T5_8+a9*T5_9+a10*T5_10
             T5_1 | b1*t1 ; T5_2 | b2*t1 ; T5_3 | b3*t1 ; T5_4 | b4*t1
             T5_5 | b5*t1 ; T5_6 | b6*t1 ; T5_7 | b7*t1 ; T5_8 | b8*t1
             T5_9 | b9*t1 ; T5_10 | b10*t1
             F1 ~~ 1*F1
             #**** T2
             F2=~ a1*T6_1+a2*T6_2+a3*T6_3+a4*T6_4+a5*T6_5+a6*T6_6+a7*T6_7+a8*T6_8+a9*T6_9+a10*T6_10
             T6_1 | b1*t1 ; T6_2 | b2*t1 ; T6_3 | b3*t1 ; T6_4 | b4*t1
             T6_5 | b5*t1 ; T6_6 | b6*t1 ; T6_7 | b7*t1 ; T6_8 | b8*t1
             T6_9 | b9*t1 ; T6_10 | b10*t1
             F2 ~~ NA*F2
             F2 ~ 1
             #*** covariance
             F1 ~~ F2
                "

mod3lav <- lavaan::cfa( data=dat, model=lavmodel,
                        std.lv=TRUE, ordered=colnames(dat), parameterization="theta")
summary(mod3lav, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)'

#*************************************************
# Model 4: Linking with items of different item slope groups
#*************************************************

# dataset for T1
dat1 <- dat[, grep( "T5", colnames(dat) ) ]
colnames(dat1) <- gsub("T5","", colnames(dat1) )
# dataset for T2
dat2 <- dat[, grep( "T7", colnames(dat) ) ]
colnames(dat2) <- gsub("T7","", colnames(dat2) )

# 2PL model with slope groups T1

mod1 <- sirt::rasch.mml2(dat1, est.a = rep(1, ncol(dat1)))
summary(mod1)

# 2PL model with slope groups T2
mod2 <- sirt::rasch.mml2( dat2, est.a = rep(1, ncol(dat1)))
summary(mod2)

#------- Link 1: Haberman Linking
# collect item parameters
dfr1 <- data.frame( "study1", mod1$item$item, mod1$item$a, mod1$item$b )
dfr2 <- data.frame( "study2", mod2$item$item, mod2$item$a, mod2$item$b )
colnames(dfr2) <- colnames(dfr1) <- c("study", "item", "a", "b" )
itempars <- rbind( dfr1, dfr2 )
# Linking
link1 <- sirt::linking.haberman(itempars=itempars)
summary(link1)
# R2 value close to 1 indicates strong invariance btw time points
# sqrtU2 values are root mean sq differences and should be small for invariance
# based on result: high R2 and low sqrtU2 and small DIF = high invariance, none excluded

#------- Link 2: Invariance alignment method
# create objects for invariance.alignment
nu <- rbind(c(mod1$item$thresh), c(mod2$item$thresh))
lambda <- rbind(c(mod1$item$a), c(mod2$item$a))
colnames(lambda) <- colnames(nu) <- c('s1','s2','s3','s4','s5','s6','s7','s8','s9','s10')
rownames(lambda) <- rownames(nu) <- c("T1", "T2")
# Linking
link2a <- sirt::invariance.alignment( lambda, nu )
summary(link2a)

