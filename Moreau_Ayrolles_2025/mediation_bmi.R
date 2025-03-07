library(tidySEM)
library(lavaan)
library(dplyr)
library(combat.enigma)

# Load file
df_global <- read.table(paste0('file.tsv'), sep = '\t', header = TRUE)
df_sub <- subset(df_global,!is.na(BMI_Z))

# Select only TD and cases of interest
case = c('case')
DATA <- subset(df_sub, principal_diag %in% c("TD", case))

# Apply Combat before the mediation

# Z score based on TD
metrics <- tail(colnames(DATA), 5)
for (m in metrics)  {
  only_td <- DATA[which(DATA$Diag == "TD"),]
  DATA[,paste0(m, "_scaled")] <- (DATA[,m]-mean(only_td[,m]))/sd(only_td[,m])
}

# Mediation
metrics_for_mediation <- tail(colnames(DATA), 5)
for (metric in metrics_for_mediation) {
  set.seed(24) 
  
  mediation_model <- paste0('
    # Regressions
    BMI_Z  ~ a * Diag + age_at_scan + Sex                          
    ', metric,' ~ c * Diag + b * BMI_Z  + age_at_scan + Sex
  
    # Indirect effect (a * b)
    indirect := a * b
  
    
    # Total effect (c + indirect)
    total := c + (a * b)
  ')
  
  mediation_results <- sem(mediation_model, data = DATA, se='boot', bootstrap=500)
  med <- summary(mediation_results, standardized = TRUE, fit.measures = TRUE, ci=TRUE)
  
  # Extract effects
  effect_direct <- med$pe$est[4]
  effect_indirect <-  med$pe$est[16]
  effect_total <-  med$pe$est[17]
  
  # Compute percentages
  percent_direct <- (abs(effect_direct) / (abs(effect_direct) + abs(effect_indirect))) * 100
  percent_indirect <- (abs(effect_indirect) / (abs(effect_direct) + abs(effect_indirect))) * 100
  percent_total <- percent_direct + percent_indirect

}