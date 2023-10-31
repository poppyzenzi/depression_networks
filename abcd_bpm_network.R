library(dplyr)
#library(psychTools)
library(qgraph)
library(tidyr)
library(bootnet)


# node strength = how well directly connected to others
# node closeness = how well indirectly connected to other nodes
# betweenness = how important in average path between two others

## pairwise Markov random field (PMRF)
# binary data = Ising model, multivariate normal density = GGM, 
# non-nomral conitnuous = requires nonparanomral transformation before estimating GGM)
# mixed graphical models can estimate PMRF w both cont and cat vars 

'9 I feel worthless or inferior
11 I am too fearful or anxious	
12 I feel too guilty	
13 I am self-conscious or easily embarrassed	
18 I am unhappy, sad, or depressed	
19 I worry a lot'

labels <- c("worthless","anxious","guilty","self-conscious","unhappy","worry")

bpm <- read.csv('/Volumes/igmm/GenScotDepression/data/abcd/release5.0/core/mental-health/mh_y_bpm.csv')
symptoms = c('bpm_9_y','bpm_11_y','bpm_12_y','bpm_13_y','bpm_18_y','bpm_19_y')
int_cols = c('src_subject_id','eventname','bpm_y_scr_internal_r', symptoms)

# select cols for internalising, only for values 0-2 in symptoms (valid)
bpm_int <- bpm %>% 
  subset(select = int_cols) %>%
  filter(if_all(all_of(symptoms), ~ . >=0 &if_all(all_of(symptoms), ~ . <=2))) %>%
  mutate(`eventname` = dplyr::recode(`eventname`,
                              `baseline_year_1_arm_1`="0",
                              `6_month_follow_up_arm_1`="1",
                              `1_year_follow_up_y_arm_1`="2",
                              `18_month_follow_up_arm_1`="3",
                              `2_year_follow_up_y_arm_1`="4",
                              `30_month_follow_up_arm_1`="5",
                              `3_year_follow_up_y_arm_1`="6",
                              `42_month_follow_up_arm_1`="7",
                              `4_year_follow_up_y_arm_1`="8"))

# make sum score
#bpm_int$sumscore <- rowSums(bpm_int[symptoms])

# save wide df with sum score only
bpm_sum_wide <- bpm_int %>%
  subset(select = c('src_subject_id','eventname','sumscore')) %>%
  pivot_wider(
  data = .,
  id_cols = src_subject_id,  # Identifier column
  names_from = eventname,    # Column to spread
  values_from = c(sumscore)  # Values to fill the columns
)

colnames(bpm_sum_wide) <- c('src_subject_id','t1','t2','t3','t4','t5','t6','t7','t8') # rename cols
bpm_sum_wide[is.na(bpm_sum_wide)] <- -9999 # missing

# now save
write.table(bpm_sum_wide, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/bpm_sum_wide.txt") # now save 

#### symptom frequency table #####
freq_dat <- bpm_int
colnames(freq_dat) <- c("id","time", "sum", unlist(labels))

frequency_table <- freq_dat %>%
  gather(symptom, value, `worthless`:`worry`) %>% # gather symptoms into key value pairs 
  group_by(time, symptom) %>% # group by time/sweep
  summarize(count = sum(value)) %>%
  spread(time, count, fill = 0)

print(frequency_table)
setwd("/Volumes/igmm/GenScotDepression/users/poppy/abcd/symptom_data")
write.csv(frequency_table, file="abcd_symptom_frequencies.csv")

##########################

# choose time point
timepoint = bpm_int$eventname == 3

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
differenceTest(boot1,"bpm_19_y","bpm_13_y","strength")



####################################
############# for mplus  ###########
####################################

# longtiudinal measurement invariance format (LMI)

bpm_wide <- pivot_wider(
  data = bpm_int,
  id_cols = src_subject_id,  # Identifier column
  names_from = eventname,    # Column to spread
  values_from = c(bpm_9_y, bpm_11_y, bpm_12_y, bpm_13_y, bpm_18_y, bpm_19_y)  # Values to fill the columns
)

bpm_wide[is.na(bpm_wide)] <- -9999

write.table(bpm_wide, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/
            Edinburgh/gmm/gmm_abcd/mplus_data/bpm_symptoms_wide.dat", col.names = FALSE)

bpm_wide_3point <- bpm_wide[,c(2:49)]

write.table(bpm_wide_3point, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/bpm_symptoms_wide_3point.txt", col.names = FALSE)


# for RI-CLPM network (make symptoms binary)
bpm_binary <- bpm_wide %>%
  mutate(across(2:49, ~case_when(
    . == 0 ~ 0,
    . == 1 ~ 1,
    . == 2 ~ 1,
    TRUE ~ .
  )))

# checks - should be 11726 (including all missing)
length(unique(bpm_binary$src_subject_id))

# drop all irrelevant cols
bpm_binary_symptoms_only <- bpm_binary[,c(2:49)]

write.table(bpm_binary_symptoms_only, "/Users/poppygrimes/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Edinburgh/networks/mplus_ri_clpm/bpm_symptoms_binary.dat", col.names = FALSE, row.names = FALSE)






