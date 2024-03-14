## MDD risk subgroups
library(dplyr)

# load in dat, select IID and PRS only
setwd('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/projects/prs/prs_alspac_OUT/all_thresholds')
mdd <- read.table('alspac_mdd_ipsych_prs_0816.t2.best',header = TRUE)
mdd <- mdd[,c(2,4)]

# Calculate quintiles
quintiles <- quantile(mdd$PRS, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1))

# Define labels for categories
labels <- c("very low", "low", "middle", "high", "very high")

# Create a new column with category labels
mdd$category <- cut(mdd$PRS, breaks = quintiles, labels = labels, include.lowest = TRUE)

# Print the first few rows to verify
head(mdd)

write.table(mdd, '/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/projects/networks/depression_networks/network_temperature/alspac_quintiles.txt')
