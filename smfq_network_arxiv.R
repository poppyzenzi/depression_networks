library(bootnet)
library(dplyr)
library(qgraph)
library(haven)

# node strength = how well directly connected to others
# node closeness = how well indirectly connected to other nodes
# betweenness = how important in average path between two others

## pairwise Markov random field (PMRF)
# binary data = Ising model, multivariate normal density = GGM, 
# non-nomral conitnuous = requires nonparanomral transformation before estimating GGM)
# mixed graphical models can estimate PMRF w both cont and cat vars 


# require smfq in long format [id, timepoint, sum score]
smfq <- read_dta('/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/gmm/gmm_alspac/alspac_dep_long.dat')

# checks
length(unique(smfq$id)) # should be 15,645

# select first 4 time points
smfq <- filter(smfq, smfq$time < 5)

# choose time point
timepoint = (smfq$time == 1)


# run bootnet for timepoint and symptom cols 
results <- bootnet(bpm_int[timepoint, 4:9], 
                   nBoots = 100, default = "ggmModSelect")

# plot network
plot(results$sample)

# ===============================


## steps for assessing accuracy of psychological network structures (Fried tutorial paper)
# a) estimate CIs on edge weights. Use bootstrapping to repeatedly estimate model under sampled data and estimate statistic of interest.
# b) assess stability of centrality indices in subsets of cases
# c) test for sig diff btw edge weights and centrality indices

# symptoms are ordinal: compute polychoric correlation matrix as input using qgraph::cor_auto


# estimate network
Network <- estimateNetwork(bpm_int[timepoint, 4:9],default="EBICglasso")

# plot estimated network
plot(Network, layout='spring', labels=TRUE)
# plot centrality indices
centralityPlot(Network, include = c("Betweenness","Closeness","Strength"))


## non-parametric bootstrap (important to use for ordinal data)
boot1<-bootnet(Network, nBoots=1000, nCores=1)

# plot bootstrapped edgeCIs: indicates care to be taken interpreting edges in network 
plot(boot1, labels=FALSE, order="sample")

# plot significant differences (alpha=0.05) of edges: 
plot(boot1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
# Plot significant differences (alpha=0.05) of nodestrength: 
plot(boot1,"strength",plot="difference")

# test for difference in strength between node x and y
x <-
y <- 
differenceTest(boot1,x,y,"strength")
