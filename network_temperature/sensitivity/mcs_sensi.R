


setwd('/Volumes/igmm/GenScotDepression/users/poppy/mcs')
mcs_symptoms <-read.table('sdq_symp_4waves_beforeQC.txt')
colnames(mcs_symptoms) <- c("ID", 
                  "s4PSDHS00", "s4PSDMW00", "s4PSDUD00", "s4PSDNC00", "s4PSDFE00",
                  "s4PSDSP00", "s4PSDGF00","s4PSDLC00", "s4PSDPB00","s4PSDGB00",
                  
                  "s5PSDHS00", "s5PSDMW00", "s5PSDUD00", "s5PSDNC00", "s5PSDFE00",
                  "s5PSDSP00", "s5PSDGF00","s5PSDLC00", "s5PSDPB00","s5PSDGB00",
                  
                  "s6PSDHS00","s6PSDMW00", "s6PSDUD00", "s6PSDNC00", "s6PSDFE00",
                  "s6PSDSP00","s6PSDGF00", "s6PSDLC00", "s6PSDPB00", "s6PSDGB00",
                  
                  "s7PSDHS00", "s7PSDMW00", "s7PSDUD00","s7PSDNC00", "s7PSDFE00",
                  "s7PSDSP00", "s7PSDGF00", "s7PSDLC00", "s7PSDPB00", "s7PSDGB00",
                  
                  "sex"
)
## QC SYMPTOM SCORES
# -1 = NA, 1 = not true, 2 = somewhat true, 3 = certainly true, 4 = blank
# symptoms scoring -1 (NA)
# males 1, females 2 --> males 0, females 1
mcs_symptoms_QC <- mcs_symptoms %>%
  replace(. == -1, NA) %>%
  replace(. == 4, NA) %>%
  replace(. == 1, 1) %>%
  replace(. == 2, 2) %>%
  replace(. == 3, 3) %>%
  replace(. == -9, NA) %>% # refusal
  replace(. == -8, NA) # don't know


# first remove labels to prevent conflicting values error
mcs_symptoms_QC <- zap_labels(mcs_symptoms_QC)

# make long
mcs_symptoms_QC_long <- mcs_symptoms_QC %>%
  pivot_longer(
    cols = starts_with("s4") | starts_with("s5") | starts_with("s6") | starts_with("s7"),
    names_to = c("Sweep", ".value"),
    names_pattern = "s(\\d+)([A-Z]+.*)"
  )

labels <- c("malaise", "worries", "unhappy", "anxiety", "fears", "solitary", "friends*", "liked*","bullied","adult-oriented")
# * reverse scored
colnames(mcs_symptoms_QC_long) <- c("id","sex","sweep", unlist(labels))

## filter for only valid symptom scores
mcs_symptoms_QC_long <- mcs_symptoms_QC_long %>%
  filter(if_all(.cols = all_of(names(mcs_symptoms_QC_long)[4:ncol(mcs_symptoms_QC_long)]), ~ . >= 0 & . <= 3))

data_subset <- as.data.frame(mcs_symptoms_QC_long[, 4:13])
par(mfrow=c(3,4)) 

for (i in 1:ncol(data_subset)) {
  hist(data_subset[, i], breaks = 0:3, main = colnames(data_subset)[i], xlab = "Values", cex.main=0.9)
}

main_title <- "Histograms of SDQ Scores in MCS"
mtext(main_title, outer = TRUE, cex = 0.8, font = 2, line = -1)

par(mfrow=c(1,1)) 
