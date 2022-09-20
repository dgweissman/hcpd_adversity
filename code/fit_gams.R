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
require(visreg)

#######################
## Read in HCPD data ##
#######################
# hcpd_data <- read.csv("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/subject_demographics/n925_hcpd_myelin_May2022_dwcleaned.csv")
hcpd_data <- read.csv("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/subject_demographics/n925_hcpd_myelin_AS.csv")

## Find subject with errant p_grade score
error_index <- which(hcpd_data$p_grade_cat==77)
hcpd_data$subject_id[error_index]
hcpd_data$p_grade_cat[error_index] <- NA
hcpd_data$p_grade[error_index] <- NA
hcpd_data$p_grade_cat <- droplevels.factor(hcpd_data$p_grade_cat)
table(hcpd_data$p_grade_cat)

## Set variables
hcpd_data$sex<- as.factor(hcpd_data$sex)
hcpd_data$site <- as.factor(hcpd_data$site)
hcpd_data$SES_RLVL <- as.factor(hcpd_data$SES_RLVL)

## Remove missing data
na_idx <- which(is.na(hcpd_data$p_grade_cat))
length(na_idx) 

ses_df <-  hcpd_data[-na_idx, ]

## Index columns with regional myelin measurements
nreg <- 360
myelin_col_index <- grep("myelin_glasser_v", colnames(hcpd_data))
## Index columns with network myelin measurements
net_col_index <- grep("_myelin", colnames(hcpd_data))

## Regional labels from Glasser parcellation
glasser_parcel_labels <- read.table("/ncf/hcp/data/analyses/myelin/parcellations/glasser360NodeNames.txt",  header=FALSE)

## Test GAM
test_gam <- gam(mean_wholebrain_T1wT2w ~ p_grade + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=hcpd_data, method="REML")
summary(test_gam)$p.table

################################################
## Estimate p_grade effect on regional myelin ##
################################################
glasser_gam_results <- NULL
mod <- NULL

## Specify covariates
full_covars=" ~  p_grade + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(ses_df[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(ses_df[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
glasser_gam_models <- lapply(m, function(x) {gam(formula = x, data=ses_df, method="REML")})
glasser_reduced_gam_models <- lapply(m2, function(x) {gam(formula = x, data=ses_df, method="REML")})
glasser_age_resid <- lapply(glasser_reduced_gam_models, function(x) { resid(x)})

glasser_p_grade_pvals <- lapply(glasser_gam_models, function(x) { summary(x)$p.table[2,4]})
glasser_p_grade_pvals <- unlist(glasser_p_grade_pvals)
glasser_p_grade_tstat <-unlist(lapply(glasser_gam_models, function(x) { summary(x)$p.table[2,3]}))

## Multiple Comparison correction
glasser_p_grade_pvals_corrected <- p.adjust(glasser_p_grade_pvals, method="holm")

corrected_sig_index <- which(glasser_p_grade_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(glasser_p_grade_pvals < 0.005)
length(uncorrected_sig_index)
glasser_parcel_labels[uncorrected_sig_index,1]
glasser_p_grade_pvals[uncorrected_sig_index]
glasser_p_grade_tstat[uncorrected_sig_index]

hist(glasser_p_grade_tstat, col="slategray3")

write.table(glasser_p_grade_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_pgrade_tstat.txt", row.names=FALSE, col.names=FALSE)

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
glasser_reduced_gam_models <- lapply(m2, function(x) {gam(formula = x, data=hcpd_data, method="REML")})
glasser_age_resid <- lapply(glasser_reduced_gam_models, function(x) { resid(x)})

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

hist(network_logINR_tstat, col="slategray3")


# write.table(network_logINR_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_logINR_tstat.txt", row.names=FALSE, col.names=FALSE)

####################################################

#################################################
## Age*logINR Interaction Across Brain Regions ##
#################################################
test_lm <- lm(mean_wholebrain_T1wT2w ~  age*logINR + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df)
summary(test_lm)$coeff
summary(test_lm)$coeff[15,4]

## Specify covariates
m <- NULL
glasser_lm_int_models <- NULL

full_covars=" ~  age*logINR + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  logINR + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(ses_df[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(ses_df[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
glasser_lm_int_models <- lapply(m, function(x) {lm(formula = x, data=ses_df)})
glasser_age_logINR_pvals <- lapply(glasser_lm_int_models, function(x) { summary(x)$coeff[15,4]})
glasser_age_logINR_pvals <- unlist(glasser_age_logINR_pvals)
glasser_age_logINR_tstat <-unlist(lapply(glasser_lm_int_models, function(x) { summary(x)$coeff[15,3]}))

max(abs(glasser_age_logINR_tstat))
hist(glasser_age_logINR_tstat, col="slategray3")

## Multiple Comparison correction
glasser_age_logINR_pvals_corrected <- p.adjust(glasser_age_logINR_pvals, method="holm")

corrected_sig_index <- which(glasser_age_logINR_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(glasser_age_logINR_pvals < 0.005)
length(uncorrected_sig_index)
glasser_parcel_labels[uncorrected_sig_index,1]
glasser_age_logINR_pvals[uncorrected_sig_index]
glasser_age_logINR_tstat[uncorrected_sig_index]

hist(glasser_age_logINR_tstat, col="slategray3")

write.table(glasser_age_logINR_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_age_by_logINR_tstat.txt", row.names=FALSE, col.names=FALSE)

## Plot Interaction effect in AVI
v111_int_plot <- ggplot(data=logINR_df, aes(x=age, y=myelin_glasser_v111)) +
  geom_smooth(data=logINR_df, aes(x=age, y=myelin_glasser_v111, col=logINR_cat2), method="lm")
v111_int_plot  + gbtheme

###################################################
## Age*logINR Interaction Across Brain Networks ##
###################################################
test_lm <- lm(Frontoparietal_myelin ~  age*logINR + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df)
summary(test_lm)$coeff
summary(test_lm)$coeff[15,4]

## Specify covariates
m <- NULL
full_covars=" ~  age*logINR + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  logINR + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(ses_df[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(ses_df[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
network_lm_int_models <- lapply(m, function(x) {lm(formula = x, data=ses_df)})
network_age_logINR_pvals <- lapply(network_lm_int_models, function(x) { summary(x)$coeff[15,4]})
network_age_logINR_pvals <- unlist(network_age_logINR_pvals)
network_age_logINR_tstat <-unlist(lapply(network_lm_int_models, function(x) { summary(x)$coeff[15,3]}))

max(abs(network_age_logINR_tstat))
hist(network_age_logINR_tstat, col="slategray3")

## Multiple Comparison correction
network_age_logINR_pvals_corrected <- p.adjust(network_age_logINR_pvals, method="holm")

corrected_sig_index <- which(network_age_logINR_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(network_age_logINR_pvals < 0.005)
length(uncorrected_sig_index)
network_age_logINR_pvals[uncorrected_sig_index]
network_age_logINR_tstat[uncorrected_sig_index]

#########################################
## Test GAM using mean cortical myelin ##
#########################################
require(gratia)

test_gam <- gam(mean_wholebrain_T1wT2w ~  s(age, by=p_grade_cat) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor,data=ses_df)

## Visual examination of spline interaction model
gratia::draw(test_gam)

## Evaluate differences of factor smooth interaction
factor_smooth_diff <- gratia::difference_smooths(test_gam, smooth = "s(age)", n = 100, ci_level = 0.95)
draw(factor_smooth_diff)

##################################################
## Age*p_grade Interaction Across Brain Regions ##
##################################################
test_lm <- lm(mean_wholebrain_T1wT2w ~  age*p_grade + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df)
summary(test_lm)$coeff
summary(test_lm)$coeff[15,4]

## Specify covariates
m <- NULL
full_covars=" ~  age*p_grade + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  p_grade + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
glasser_lm_int_models <- lapply(m, function(x) {lm(formula = x, data=ses_df)})
glasser_age_pgrade_pvals <- lapply(glasser_lm_int_models, function(x) { summary(x)$coeff[15,4]})
glasser_age_pgrade_pvals <- unlist(glasser_age_pgrade_pvals)
glasser_age_pgrade_tstat <-unlist(lapply(glasser_lm_int_models, function(x) { summary(x)$coeff[15,3]}))

max(abs(glasser_age_pgrade_tstat))
hist(glasser_age_pgrade_tstat, col="slategray3")

which(glasser_age_pgrade_tstat <  -2.7)

## Multiple Comparison correction
glasser_age_pgrade_pvals_corrected <- p.adjust(glasser_age_pgrade_pvals, method="holm")

corrected_sig_index <- which(glasser_age_pgrade_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(glasser_age_pgrade_pvals < 0.005)
length(uncorrected_sig_index)
glasser_parcel_labels[uncorrected_sig_index,1]
glasser_age_pgrade_pvals[uncorrected_sig_index]
glasser_age_pgrade_tstat[uncorrected_sig_index]

write.table(glasser_age_pgrade_tstat, "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_age_by_pgrade_tstat.txt", row.names=FALSE, col.names=FALSE)

###################################
## VISUALIZE INTERACTION EFFECTS ##
###################################


## ggplot theme (aesthetics) 
apal <- paste0('#', c('000000', 'EAE3D8', 'FFFFFF', 'FFD378', '424C6D'))

gbtheme <- theme_bw() +  
  theme(text = element_text(size = 20),
        panel.background = element_rect(fill = apal[[3]], size = 0, color = apal[[2]]),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = apal[[2]], size = 0),
        strip.text = element_text(color = '#222222'),
        axis.text =  element_text(color = apal[[1]]), axis.title = element_text(color = apal[[1]]),
        axis.ticks.length=unit(0.25, "cm"))

## Create 2-level factor for p_grade (high/low)
quantile(ses_df$p_grade, probs = c(0.33, 0.66))

length(which(ses_df$p_grade > 18 & ses_df$p_grade < 20))
length(which(ses_df$p_grade > 18))

low_idx <- which(ses_df$p_grade < 18)
high_idx <- which(ses_df$p_grade >= 18)

## Create new variable for interaction visualization
ses_df$p_grade_cat2 <- 0
ses_df$p_grade_cat2[low_idx] <- 1  
ses_df$p_grade_cat2[high_idx] <- 2

ses_df$p_grade_cat2 <- as.factor(ses_df$p_grade_cat2)

## Right Anterior ventral insula
v111_lm <- lm(myelin_glasser_v111 ~  age*p_grade + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df)

## Interaction with 2 levels of p_grade
v111_int_plot <- ggplot(data=ses_df, aes(x=age, y=myelin_glasser_v111)) +
  geom_smooth(data=ses_df, aes(x=age, y=myelin_glasser_v111, col=p_grade_cat2), method="lm")
v111_int_plot  + gbtheme

## Interaction with 5 levels of p_grade
v111_int_plot2 <- ggplot(data=ses_df, aes(x=age, y=myelin_glasser_v111)) +
  geom_smooth(data=ses_df, aes(x=age, y=myelin_glasser_v111, col=p_grade_cat), method="lm")
v111_int_plot2  + gbtheme

###################################################
## Age*p_grade Interaction Across Brain Networks ##
###################################################
test_lm <- lm(Frontoparietal_myelin ~  age*p_grade + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df)
summary(test_lm)$coeff

## Specify covariates
m <- NULL
full_covars=" ~  age*p_grade + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  p_grade + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
network_lm_int_models <- lapply(m, function(x) {lm(formula = x, data=ses_df)})
network_age_pgrade_pvals <- lapply(network_lm_int_models, function(x) { summary(x)$coeff[15,4]})
network_age_pgrade_pvals <- unlist(network_age_pgrade_pvals)
network_age_pgrade_tstat <-unlist(lapply(network_lm_int_models, function(x) { summary(x)$coeff[15,3]}))

max(abs(network_age_pgrade_tstat))

## Multiple Comparison correction
network_age_pgrade_pvals_corrected <- p.adjust(network_age_pgrade_pvals, method="holm")

corrected_sig_index <- which(network_age_pgrade_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(network_age_pgrade_pvals < 0.005)
length(uncorrected_sig_index)
network_age_pgrade_pvals[uncorrected_sig_index]
network_age_pgrade_tstat[uncorrected_sig_index]

hist(network_age_pgrade_tstat, col="slategray3")

##############################################
## Age*SES Interaction Across Brain Regions ##
##############################################

test_lm <- lm(mean_wholebrain_T1wT2w ~  age*SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=ses_df_2)
summary(test_lm)$coeff

## Specify covariates
m <- NULL
full_covars=" ~  age*SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(ses_df_2[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(ses_df_2[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
glasser_lm_int_models <- lapply(m, function(x) {lm(formula = x, data=ses_df_2)})
glasser_age_lowSES_pvals <- lapply(glasser_lm_int_models, function(x) { summary(x)$coeff[16,4]})
glasser_age_lowSES_pvals <- unlist(glasser_age_lowSES_pvals)
glasser_age_lowSES_tstat <-unlist(lapply(glasser_lm_int_models, function(x) { summary(x)$coeff[16,3]}))

max(abs(glasser_age_lowSES_tstat))
hist(glasser_age_lowSES_tstat, col="slategray3")

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

## Create SES df (Drop 62 "UKNOWN" and "NA" cases)
ses_df_2 <- subset(hcpd_data, hcpd_data$SES_RLVL=="HIGH" | hcpd_data$SES_RLVL=="MIDDLE" | hcpd_data$SES_RLVL=="LOW")
ses_df_2$SES_RLVL <- as.factor(ses_df_2$SES_RLVL)
levels(ses_df_2$SES_RLVL)
ses_df_2$SES_RLVL <- droplevels(ses_df_2$SES_RLVL)
levels(ses_df_2$SES_RLVL)


## Specify covariates
m <- NULL
full_covars=" ~  age*SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  SES_RLVL + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(ses_df_2[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(ses_df_2[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
network_lm_int_models <- lapply(m, function(x) {lm(formula = x, data=ses_df_2)})
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
R_47m_plot <- ggplot(data=ses_df_2, aes(x=age, y=myelin_glasser_v66)) +
  geom_smooth(data=ses_df_2, aes(x=age, y=myelin_glasser_v66, col=SES_RLVL), method="gam", formula=y~ s(x))
R_47m_plot 
R_47m_plot + gbtheme

## Right_A4 interaction
R_A4_plot <- ggplot(data=ses_df_2, aes(x=age, y=myelin_glasser_v175)) +
  geom_smooth(data=ses_df_2, aes(x=age, y=myelin_glasser_v175, col=SES_RLVL), method="gam", formula=y~ s(x))
R_A4_plot
R_A4_plot + gbtheme