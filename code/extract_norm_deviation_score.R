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

# Create empty data-frame for normative difference scores
regional_norm_diff_Zscores <- data.frame(matrix(nrow = nsub, ncol = length(response_vars)))

# Load gaussian process model outputs ----
regional_z_scores_list <- lapply(1:length(response_vars), function(i) {

  ## Myelin
  right_side <- '~ gp(age, k=10, c=5/4) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'
  left_side <- response_vars[[i]]
  
  ## Thickness
  # right_side <- '~ gp(Age, k=10, c=5/4) + Sex + Scanner'
  # left_side <- response_vars[[i]]
  
  mf_form <- as.formula(paste0(left_side, gsub('gp\\((\\w+).*\\)', '\\1', right_side)))
  form <- as.formula(paste0(left_side, right_side))
  
  d <- model.frame(formula = mf_form, data = hcpd_data)
  
  ## Return subject IDs for later merging
  d_row_index <- as.numeric(row.names(d))
  d$subject_id <- hcpd_data$subject_id[d_row_index]
  
  # Identify response variable 
  col_idx <- grep(response_vars[[i]], colnames(d))
  
  # Scale and center response variable 
  d[, col_idx] <- scale(d[, col_idx], center = TRUE, scale = TRUE)
  
  # Gaussian process model fit
  model_name <- paste0('gp_fit_', left_side)
  
  if(!file.exists(paste0(model_name, '.rds'))){
    gp_fit <- NULL
  } else {
    gp_fit <- readRDS(paste0(model_name, '.rds'))
  }
  
  # 10-fold cross validation ----
  if(!file.exists(paste0('kfold_gp_fit_', left_side, '.rds'))){
    future::plan(future::multiprocess, workers = 10)
    kfold1 <- brms::kfold(gp_fit, K = 10, save_fits = TRUE, chains = 4, nug = 1e-07)
    saveRDS(kfold1, paste0('kfold_gp_fit_', left_side, '.rds'))
  } else {
    kfold1 <- readRDS(paste0('kfold_gp_fit_', left_side, '.rds'))
  }
  
  ## Function for calculating out-of-sample normative stats 
  oos_norm_Z <- function(pe_sample, sigma_nj) {
    qstats <- quantile(pe_sample, probs = c(.025, .5, .975))
    m <- mean(pe_sample) #y_obs - y_hat
    se <- sd(pe_sample) #\sigma_ij
    Z <- m / sqrt( se^2 + sigma_nj^2 )
    stats <- c(qstats, mean = m, se = se, Z = Z)
    return(stats)
  }
  
  ## For each fold in kfold output, calculate normative scores for "left-out" subjects 
  norm_diffs <- apply(kfold1$fits[, 1:2], 1, function(arow, data = kfold1$data, d_sid = d){
    afit <- arow[['fit']]
    omitted <- arow[['omitted']]
    nd <- data[omitted, ]
    nd$id <- omitted
    nd$subject_id <- d_sid$subject_id[omitted]
    pe <- brms::predictive_error(afit, newdata = nd, nug = 1e-07)
    sigma <- brms::posterior_summary(afit, pars = 'sigma')[,'Estimate']
    stats <- apply(pe, 2, function(pe_col) {
      s <- oos_norm_Z(pe_sample = pe_col, sigma_nj = sigma)
      return(s)
    })
    return(cbind(as.data.frame(t(stats)), nd))
  })
  
  ## Merge norm difference scores from each fold
  norm_diffs_df <- do.call(rbind, norm_diffs)
  
  ## Subset response variable and subject_id for merging. We just want the
  ## subject ID and the Z score.
  vars <- c("Z", "subject_id")
  slim_df <- norm_diffs_df[vars]
  #We also want to add a column id that identifies the region.
  slim_df$region_id <- gsub('myelin_glasser_', '', left_side)
  return(slim_df)
})

#rbind all the data into a single data table. I prefer data.table for combining
#alot of data like this as it is always faster.
system.time({regional_z_scores_dt <- data.table::rbindlist(regional_z_scores_list)})

system.time({regional_z_scores_dt_wide <- data.table::dcast(regional_z_scores_dt, ... ~ region_id, value.var = 'Z')})

#how big are these objects?
format(object.size(regional_z_scores_dt), units = 'MB')
format(object.size(regional_z_scores_dt_wide), units = 'MB')

#The wide version is smaller on disk.
write.csv(regional_z_scores_dt_wide, '/ncf/hcp/data/analyses/myelin/Adversity_Project/results/normative_modeling/regional/zscore/n925_regional_myelin_normative_Zscores_wide.csv', row.names=FALSE)
# write.csv(regional_z_scores_dt, '/ncf/hcp/data/analyses/myelin/normative_modeling_sandbox/Dec2021/myelin/regional_Zscores/n628_regional_myelin_normative_Zscores.csv')

#######################
### READ IN ZSCORES ###
#######################
regional_z_scores_dt_wide <- read.csv("/ncf/hcp/data/analyses/myelin/Adversity_Project/results/normative_modeling/regional/zscore/n925_regional_myelin_normative_Zscores_wide.csv")

## Change names ##
##################
nreg <- 360
new_df <- regional_z_scores_dt_wide[,3:362]
subject_id <- regional_z_scores_dt_wide$subject_id
for(i in 1:nreg) {
  tmp_name <- names(new_df)[i]
  names(new_df)[i] <- paste('myelin_glasser_', tmp_name,'_normDiff_zscore', sep="")
  }

new_df <- cbind(subject_id, new_df)
norm_df <- merge(hcpd_data, new_df, by="subject_id")

## export df
write.csv(norm_df, "/ncf/hcp/data/analyses/myelin/Adversity_Project/data/n925_hcpd_normDiff_myelin.csv")

##############################

# Identify response variables
col_idx <- grep("normDiff", colnames(norm_df))

test_lm <- lm(myelin_glasser_v100_normDiff_zscore ~  logINR + age, data=norm_df)
summary(test_lm)$coeff

## Specify covariates
m <- NULL
full_covars=" ~  logINR + age"  #+ sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   

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
