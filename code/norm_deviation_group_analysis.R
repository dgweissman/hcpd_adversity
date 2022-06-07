rm(list = ls())

## Check which modules are loaded
.libPaths(c('/ncf/hcp/data/analyses/myelin/verse-cmdstan-container', .libPaths()))
.libPaths()

# R Setup ----
library(brms)
library(data.table)
library(ggplot2)
library(readr)
library(future)
library(stringr)

## Set working directory to where model outputs are generated
setwd('/ncf/hcp/data/analyses/myelin/Adversity_Project/results/normative_modeling/regional')

## Load HCP Data
hcpd_data <- read.csv("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/subject_demographics/n925_hcpd_myelin_May2022.csv")
hcpd_data$sex<- as.factor(hcpd_data$sex)
hcpd_data$site <- as.factor(hcpd_data$site)

# Define response variables
response_vars <- c(paste0('myelin_glasser_v', 1:360)) 

# Define number of subjects
nsub <- dim(hcpd_data)[1]

#################################
### READ IN NORMATIVE ZSCORES ###
#################################
regional_z_scores_dt_wide <- read.csv("/ncf/hcp/data/analyses/myelin/Adversity_Project/results/normative_modeling/regional/zscore/n925_regional_myelin_normative_Zscores_wide.csv")

##################
## Change names ##
##################
nreg <- 360
new_df <- regional_z_scores_dt_wide[,2:361]
subject_id <- regional_z_scores_dt_wide$subject_id

for(i in 1:nreg) {
  tmp_name <- names(new_df)[i]
  names(new_df)[i] <- paste('myelin_glasser_', tmp_name,'_normDiff_zscore', sep="")
}

names(new_df)

new_df <- cbind(subject_id, new_df)
norm_df <- merge(hcpd_data, new_df, by="subject_id")

## export df
write.csv(norm_df, "/ncf/hcp/data/analyses/myelin/Adversity_Project/data/n925_hcpd_normDiff_myelin.csv", row.names=FALSE)

############################################################


test_lm <- lm(myelin_glasser_v100_normDiff_zscore ~  logINR + age, data=norm_df)
summary(test_lm)

## Create dataframe for linear regression output
regional_normDiff_results<- data.frame(matrix(nrow = nreg, ncol = 3))
colnames(regional_normDiff_results) <- c("lm_betas", "lm_tstats", "lm_pvals")


for(i in 1:nreg){
  ## Run Linear Regression Model
  left_side <- paste('myelin_glasser_v', i, '_normDiff_zscore ~', sep='')
  right_side <- " logINR + age"
  tmp_lm <- lm(as.formula(paste(left_side, right_side, sep='')), data= norm_df)
  ## Output results
  regional_normDiff_results$lm_betas[i] <- summary(tmp_lm)$coeff[2,1]
  regional_normDiff_results$lm_tstats[i] <- summary(tmp_lm)$coeff[2,3]
  regional_normDiff_results$lm_pvals[i] <- summary(tmp_lm)$coeff[2,4]
  # regional_normDiff_results$adj.rsq[i] <- summary(tmp_lm)$adj.r.squared
}

hist(regional_normDiff_results$lm_tstats, col="slategray3")

normDiff_logINR_tstat <- regional_normDiff_results$lm_tstats

write.table(normDiff_logINR_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_logINR_normDiff_myelin_ageCov_tstat.txt", row.names=FALSE, col.names=FALSE)



##########################################################
# Identify response variables
col_idx <- grep("normDiff", colnames(norm_df))

## Specify covariates
m <- NULL
full_covars=" ~  logINR"  # + age + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply lm across brain regions (columns of data-frame)
m <- lapply(names(norm_df[,col_idx]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract model estimates
lm_normDiff_models <- lapply(m, function(x) {lm(formula = x, data=norm_df)})

normDiff_logINR_pvals <- unlist(lapply(lm_normDiff_models, function(x) { summary(x)$coeff[2,4]}))
normDiff_logINR_tstat <-unlist(lapply(lm_normDiff_models, function(x) { summary(x)$coeff[2,3]}))

max(abs(normDiff_logINR_tstat))

## Multiple Comparison correction
normDiff_logINR_pvals_corrected <- p.adjust(normDiff_logINR_pvals, method="holm")

corrected_sig_index <- which(normDiff_logINR_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(normDiff_logINR_pvals < 0.005)
length(uncorrected_sig_index)
glasser_parcel_labels[uncorrected_sig_index,1]
normDiff_logINR_pvals[uncorrected_sig_index]
normDiff_logINR_tstat[uncorrected_sig_index]

write.table(normDiff_logINR_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_logINR_normDiff_myelin_tstat.txt", row.names=FALSE, col.names=FALSE)
