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
hcpd_data <- read.csv("n925_hcpd_myelin_AS.csv")

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
glasser_parcel_labels <- read.table("glasser360NodeNames.txt",  header=FALSE)

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
glasser_p_grade_coefs <- unlist(lapply(glasser_gam_models, function(x) { summary(x)$p.table[2,1]}))

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

#############################################
# Estimate p_grade effect on NETWORK myelin ##
#############################################
network_gam_models <- NULL
m <- NULL

## Specify covariates
full_covars=" ~  p_grade + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
network_gam_models <- lapply(m, function(x) {gam(formula = x, data=hcpd_data, method="REML")})
network_pgrade_pvals <- lapply(network_gam_models, function(x) { summary(x)$p.table[2,4]})
network_pgrade_pvals <- unlist(network_pgrade_pvals)
network_pgrade_tstat <-unlist(lapply(network_gam_models, function(x) { summary(x)$p.table[2,3]}))

## Multiple Comparison correction
network_logINR_pvals_corrected <- p.adjust(network_logINR_pvals, method="holm")

corrected_sig_index <- which(network_logINR_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(network_logINR_pvals < 0.005)
length(uncorrected_sig_index)
# glasser_parcel_labels[uncorrected_sig_index,1]
network_logINR_pvals[uncorrected_sig_index]
network_logINR_tstat[uncorrected_sig_index]



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


### Null Inference ROPE Tests ###
# set weakly informative priors
prior1 <- c(set_prior("normal(0, 10)", class = "b"), set_prior("normal(0, 10)", class = "Intercept"))
## whole brain logINR ROPE Test ## this will take >20 minutes
## Test GAM
test_gam <- gam(mean_wholebrain_T1wT2w ~ logINR + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=hcpd_data, method="REML")
summary(test_gam)$p.table
logINR.brm <- brm(mean_wholebrain_T1wT2w ~ logINR + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor,
                   prior = prior1,
                   data = hcpd_data,
                   family = gaussian(link = "identity"),
                   warmup = 1e3, iter = 1.5e4, thin = 5,
                   chains = 4, cores = 4, seed = "123",
                   control = list(adapt_delta = 0.999,
                                  max_treedepth = 20))

# extract posterior samples of population-level effects
samples.logINR.brm <- posterior_samples(logINR.brm, "logINR")

#Calculate (sd of logINR)/(sd of mean_wholebrain_T1wT2w)
stdratio<-sd(hcpd_data$logINR,na.rm = T)/sd(hcpd_data$mean_wholebrain_T1wT2w)

# ROPES set at a standardized coefficient of +/-.01 through +/-.1 in increments of 0.01
# estimate posterior probability of coef being within the ROPE intervals
ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.01 & samples.logINR.brm$b_logINR*stdratio < .01
cat("Probability that standardized beta of logINR is > -.01 & < .01 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.02 & samples.logINR.brm$b_logINR*stdratio < .02
cat("Probability that standardized beta of logINR is > -.02 & < .02 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.03 & samples.logINR.brm$b_logINR*stdratio < .03
cat("Probability that standardized beta of logINR is > -.03 & < .03 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.04 & samples.logINR.brm$b_logINR*stdratio < .04
cat("Probability that standardized beta of logINR is > -.04 & < .04 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.05 & samples.logINR.brm$b_logINR*stdratio < .05
cat("Probability that standardized beta of logINR is > -.05 & < .05 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.06 & samples.logINR.brm$b_logINR*stdratio < .06
cat("Probability that standardized beta of logINR is > -.06 & < .06 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.07 & samples.logINR.brm$b_logINR*stdratio < .07
cat("Probability that standardized beta of logINR is > -.07 & < .07 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.08 & samples.logINR.brm$b_logINR*stdratio < .08
cat("Probability that standardized beta of logINR is > -.08 & < .08 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.09 & samples.logINR.brm$b_logINR*stdratio < .09
cat("Probability that standardized beta of logINR is > -.09 & < .09 =", mean(ROPE_test), "\n")

ROPE_test <- samples.logINR.brm$b_logINR*stdratio > -.10 & samples.logINR.brm$b_logINR*stdratio < .10
cat("Probability that standardized beta of logINR is > -.10 & < .10 =", mean(ROPE_test), "\n")

## whole brain p_grade ROPE Test ## this will take >20 minutes
## Test GAM
test_gam <- gam(mean_wholebrain_T1wT2w ~ p_grade + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor, data=hcpd_data, method="REML")
summary(test_gam)$p.table
p_grade.brm <- brm(mean_wholebrain_T1wT2w ~ p_grade + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor,
                  prior = prior1,
                  data = hcpd_data,
                  family = gaussian(link = "identity"),
                  warmup = 1e3, iter = 1.5e4, thin = 5,
                  chains = 4, cores = 4, seed = "123",
                  control = list(adapt_delta = 0.999,
                                 max_treedepth = 20))

# extract posterior samples of population-level effects
samples.p_grade.brm <- posterior_samples(p_grade.brm, "p_grade")

#Calculate (sd of p_grade)/(sd of mean_wholebrain_T1wT2w)
stdratio<-sd(hcpd_data$p_grade,na.rm = T)/sd(hcpd_data$mean_wholebrain_T1wT2w)

# ROPES set at +/- standardized coefficient of .01 through +/-.07 in increments of 0.01
# estimate posterior probability of coef being within the ROPE intervals
ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.01 & samples.p_grade.brm$b_p_grade*stdratio < .01
cat("Probability that standardized beta of p_grade is > -.01 & < .01 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.02 & samples.p_grade.brm$b_p_grade*stdratio < .02
cat("Probability that standardized beta of p_grade is > -.02 & < .02 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.03 & samples.p_grade.brm$b_p_grade*stdratio < .03
cat("Probability that standardized beta of p_grade is > -.03 & < .03 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.04 & samples.p_grade.brm$b_p_grade*stdratio < .04
cat("Probability that standardized beta of p_grade is > -.04 & < .04 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.05 & samples.p_grade.brm$b_p_grade*stdratio < .05
cat("Probability that standardized beta of p_grade is > -.05 & < .05 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.06 & samples.p_grade.brm$b_p_grade*stdratio < .06
cat("Probability that standardized beta of p_grade is > -.06 & < .06 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.07 & samples.p_grade.brm$b_p_grade*stdratio < .07
cat("Probability that standardized beta of p_grade is > -.07 & < .07 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.08 & samples.p_grade.brm$b_p_grade*stdratio < .08
cat("Probability that standardized beta of p_grade is > -.08 & < .08 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.09 & samples.p_grade.brm$b_p_grade*stdratio < .09
cat("Probability that standardized beta of p_grade is > -.09 & < .09 =", mean(ROPE_test), "\n")

ROPE_test <- samples.p_grade.brm$b_p_grade*stdratio > -.1 & samples.p_grade.brm$b_p_grade*stdratio < .1
cat("Probability that standardized beta of p_grade is > -.1 & < .1 =", mean(ROPE_test), "\n")

## network logINR ROPE Test ## this will take >20 minutes
network_gam_models <- NULL
m <- NULL
full_covars=" ~  logINR + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
m <- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
network.logINR.brm <- lapply(m,function(x){brm(formula=x,prior = prior1, data = hcpd_data, family = gaussian(link = "identity"), warmup = 1e3, iter = 1.5e4, thin = 5, chains = 4, cores = 4, seed = "123", control = list(adapt_delta = 0.999, max_treedepth = 20))})
network.logINR.ROPE<-data.frame(matrix(NA,nrow = 12,ncol = 11))
colnames(network.logINR.ROPE)[1]<-"network"
colnames(network.logINR.ROPE)[2:11]<-paste("pbeta",substr(seq(.01,.1,by=.01),2,4),sep="")
network.logINR.ROPE$network<-colnames(hcpd_data[,net_col_index])
for (a in 1:12){
  stdratio<-sd(hcpd_data$logINR,na.rm = T)/sd(hcpd_data[,net_col_index[a]])
  samples.logINR.brm <- posterior_samples(network.logINR.brm[a], "b_logINR")
  for (b in 1:10){
    network.logINR.ROPE[a,1+b]<-mean(samples.logINR.brm$b_logINR*stdratio > -.01*b & samples.logINR.brm$b_logINR*stdratio < .01*b)
  }
}
write.csv(network.logINR.ROPE,"network.logINR.ROPE.csv")

## network p_grade ROPE Test 
m <- NULL
## Specify covariates
full_covars=" ~  p_grade + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
m <- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
network.p_grade.brm <- lapply(m,function(x){brm(formula=x,prior = prior1, data = hcpd_data, family = gaussian(link = "identity"), warmup = 1e3, iter = 1.5e4, thin = 5, chains = 4, cores = 4, seed = "123", control = list(adapt_delta = 0.999, max_treedepth = 20))})
network.p_grade.ROPE<-data.frame(matrix(NA,nrow = 12,ncol = 11))
colnames(network.p_grade.ROPE)[1]<-"network"
colnames(network.p_grade.ROPE)[2:11]<-paste("pbeta",substr(seq(.01,.1,by=.01),2,4),sep="")
network.p_grade.ROPE$network<-colnames(hcpd_data[,net_col_index])
for (a in 1:12){
  stdratio<-sd(hcpd_data$p_grade,na.rm = T)/sd(hcpd_data[,net_col_index[a]])
  samples.p_grade.brm <- posterior_samples(network.p_grade.brm[a], "b_p_grade")
  for (b in 1:10){
    network.p_grade.ROPE[a,1+b]<-mean(samples.p_grade.brm$b_p_grade*stdratio > -.01*b & samples.p_grade.brm$b_p_grade*stdratio < .01*b)
  }
}
write.csv(network.p_grade.ROPE,"network.p_grade.ROPE.csv")

## SES Composite
# Calculate Composite
hcpd_data$SES_Composite<-(scale(hcpd_data$p_grade)+scale(hcpd_data$INR))/2
# Repeat Regional Myelin Analyses
glasser_gam_results <- NULL
mod <- NULL

## Specify covariates
full_covars=" ~  SES_Composite + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,myelin_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
glasser_gam_models <- lapply(m, function(x) {gam(formula = x, data=hcpd_data, method="REML")})
glasser_reduced_gam_models <- lapply(m2, function(x) {gam(formula = x, data=hcpd_data, method="REML")})
glasser_age_resid <- lapply(glasser_reduced_gam_models, function(x) { resid(x)})

glasser_SES_pvals <- lapply(glasser_gam_models, function(x) { summary(x)$p.table[2,4]})
glasser_SES_pvals <- unlist(glasser_SES_pvals)
glasser_SES_tstat <-unlist(lapply(glasser_gam_models, function(x) { summary(x)$p.table[2,3]}))

## Multiple Comparison correction
glasser_SES_pvals_corrected <- p.adjust(glasser_SES_pvals, method="holm")

corrected_sig_index <- which(glasser_SES_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(glasser_SES_pvals < 0.005)
length(uncorrected_sig_index)
glasser_parcel_labels[uncorrected_sig_index,1]
glasser_SES_pvals[uncorrected_sig_index]
glasser_SES_tstat[uncorrected_sig_index]

# Estimate effect on NETWORK myelin 
network_gam_models <- NULL
m <- NULL

## Specify covariates
full_covars=" ~  SES_Composite + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
reduced_covars=" ~  s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

## Apply GAM across brain regions (columns of data-frame)
m <- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})
m2<- lapply(names(hcpd_data[,net_col_index]), function(x) {as.formula(paste(x, full_covars, sep=""))})

## Extract GAM pvals 
network_gam_models <- lapply(m, function(x) {gam(formula = x, data=hcpd_data, method="REML")})
network_SES_pvals <- lapply(network_gam_models, function(x) { summary(x)$p.table[2,4]})
network_SES_pvals <- unlist(network_SES_pvals)
network_SES_tstat <-unlist(lapply(network_gam_models, function(x) { summary(x)$p.table[2,3]}))

## Multiple Comparison correction
network_SES_pvals_corrected <- p.adjust(network_SES_pvals, method="holm")

corrected_sig_index <- which(network_SES_pvals_corrected < 0.05)
length(corrected_sig_index)

uncorrected_sig_index <- which(network_SES_pvals < 0.005)
length(uncorrected_sig_index)
glasser_parcel_labels[uncorrected_sig_index,1]
network_SES_pvals[uncorrected_sig_index]
network_SES_tstat[uncorrected_sig_index]

hist(network_SES_tstat, col="slategray3")



