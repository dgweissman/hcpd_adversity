rm(list = ls())

##############################
# Load relevant Libraries ----
##############################
require(ggplot2)
require(stringr)
require(data.table)
require(readr)
require(mgcv)
require(gratia)
require(parallel)

#######################
## Read in HCPD data ##
#######################
hcpd_data <- read.csv("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/subject_demographics/n925_hcpd_myelin_May2022.csv")
## Set variables
hcpd_data$sex<- as.factor(hcpd_data$sex)
hcpd_data$site <- as.factor(hcpd_data$site)

## Index columns with regional myelin measurements
nreg <- 360
myelin_col_index <- grep("myelin_glasser_v", colnames(hcpd_data))
## Index columns with network myelin measurements
net_col_index <- grep("_myelin", colnames(hcpd_data))

## Regional labels from Glasser parcellation
glasser_parcel_labels <- read.table("/ncf/hcp/data/analyses/myelin/parcellations/glasser360NodeNames.txt",  header=FALSE)

##############################################
# Estimate logINR effect on regional myelin ##
##############################################
glasser_gam_results <- NULL
mod <- NULL

## Specify covariates
full_covars=" ~  logINR + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
glasser_gam_models <- lapply(m, function(x) {gam(formula = x, data=hcpd_data, method="REML")})
glasser_logINR_pvals <- lapply(glasser_gam_models, function(x) { summary(x)$p.table[2,4]})
glasser_logINR_pvals <- unlist(glasser_logINR_pvals)
glasser_logINR_tstat <-unlist(lapply(glasser_gam_models, function(x) { summary(x)$p.table[2,3]}))

## Multiple Comparison correction
glasser_logINR_pvals_corrected <- p.adjust(glasser_logINR_pvals, method="holm")

corrected_sig_index <- which(glasser_logINR_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(glasser_logINR_pvals < 0.005)
length(uncorrected_sig_index)
glasser_parcel_labels[uncorrected_sig_index,1]
glasser_logINR_pvals[uncorrected_sig_index]
glasser_logINR_tstat[uncorrected_sig_index]

write.table(glasser_logINR_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_logINR_tstat.txt", row.names=FALSE, col.names=FALSE)



#############################################
# Estimate logINR effect on NETWORK myelin ##
#############################################
network_gam_models <- NULL
m <- NULL

## Specify covariates
full_covars=" ~  logINR + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
network_gam_models <- lapply(m, function(x) {gam(formula = x, data=hcpd_data, method="REML")})
network_logINR_pvals <- lapply(network_gam_models, function(x) { summary(x)$p.table[2,4]})
network_logINR_pvals <- unlist(network_logINR_pvals)
network_logINR_tstat <-unlist(lapply(network_gam_models, function(x) { summary(x)$p.table[2,3]}))

## Multiple Comparison correction
network_logINR_pvals_corrected <- p.adjust(network_logINR_pvals, method="holm")

corrected_sig_index <- which(network_logINR_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(network_logINR_pvals < 0.005)
length(uncorrected_sig_index)
# glasser_parcel_labels[uncorrected_sig_index,1]
network_logINR_pvals[uncorrected_sig_index]
network_logINR_tstat[uncorrected_sig_index]

# write.table(network_logINR_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_logINR_tstat.txt", row.names=FALSE, col.names=FALSE)

####################################################

#########################################
## Test GAM using mean cortical myelin ##
#########################################
require(gratia)

ses_df <- subset(hcpd_data, hcpd_data$SES_RLVL!="UNKNOWN")
ses_df$SES_RLVL <- as.factor(ses_df$SES_RLVL)
levels(ses_df$SES_RLVL)

test_gam <- gam(mean_wholebrain_T1wT2w ~  s(age, by=SES_RLVL) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor,data=ses_df)

## Visual examination of spline interaction model
gratia::draw(test_gam)

## Evaluate differences of factor smooth interaction
factor_smooth_diff <- gratia::difference_smooths(test_gam, smooth = "s(age)", n = 100, ci_level = 0.95)
draw(factor_smooth_diff)


##############################################
## Age*SES Interaction Across Brain Regions ##
##############################################

test_lm <- lm(mean_wholebrain_T1wT2w ~  age*SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df)
summary(test_lm)$coeff

## Specify covariates
m <- NULL
full_covars=" ~  age*SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
glasser_lm_int_models <- lapply(m, function(x) {lm(formula = x, data=ses_df)})
glasser_age_lowSES_pvals <- lapply(glasser_lm_int_models, function(x) { summary(x)$coeff[16,4]})
glasser_age_lowSES_pvals <- unlist(glasser_age_lowSES_pvals)
glasser_age_lowSES_tstat <-unlist(lapply(glasser_lm_int_models, function(x) { summary(x)$coeff[16,3]}))

max(abs(glasser_age_lowSES_tstat))

## Multiple Comparison correction
glasser_age_lowSES_pvals_corrected <- p.adjust(glasser_age_lowSES_pvals, method="holm")

corrected_sig_index <- which(glasser_age_lowSES_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(glasser_age_lowSES_pvals < 0.005)
length(uncorrected_sig_index)
glasser_parcel_labels[uncorrected_sig_index,1]
glasser_age_lowSES_pvals[uncorrected_sig_index]
glasser_age_lowSES_tstat[uncorrected_sig_index]

write.table(glasser_age_lowSES_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_age_by_lowSES_tstat.txt", row.names=FALSE, col.names=FALSE)


###############################################
## Age*SES Interaction Across Brain Networks ##
###############################################
## Specify covariates
m <- NULL
full_covars=" ~  age*SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
network_lm_int_models <- lapply(m, function(x) {lm(formula = x, data=ses_df)})
network_age_lowSES_pvals <- lapply(network_lm_int_models, function(x) { summary(x)$coeff[16,4]})
network_age_lowSES_pvals <- unlist(network_age_lowSES_pvals)
network_age_lowSES_tstat <-unlist(lapply(network_lm_int_models, function(x) { summary(x)$coeff[16,3]}))

max(abs(network_age_lowSES_tstat))

## Multiple Comparison correction
network_age_lowSES_pvals_corrected <- p.adjust(network_age_lowSES_pvals, method="holm")

corrected_sig_index <- which(network_age_lowSES_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(network_age_lowSES_pvals < 0.005)
length(uncorrected_sig_index)
network_age_lowSES_pvals[uncorrected_sig_index]
network_age_lowSES_tstat[uncorrected_sig_index]

### VISUALIZE INTERACTION EFFECTS ####

v66_lm <- lm(myelin_glasser_v66 ~  age*SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df)
v66_gam <- gam(myelin_glasser_v66 ~  s(age, by=SES_RLVL) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df, method="REML")

## Right_47m interaction
ggplot(data=ses_df, aes(x=age, y=myelin_glasser_v66)) +
  geom_smooth(data=ses_df, aes(x=age, y=myelin_glasser_v66, col=SES_RLVL), method="gam", formula=y~ s(x))

## Right_A4 interaction
ggplot(data=ses_df, aes(x=age, y=myelin_glasser_v175)) +
  geom_smooth(data=ses_df, aes(x=age, y=myelin_glasser_v175, col=SES_RLVL), method="gam", formula=y~ s(x))
