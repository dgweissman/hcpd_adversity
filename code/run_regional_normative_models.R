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

## Task index for parallelizing regional models
task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))


# for (i in 1:360){

  region_id <- task_id

  set.seed(123)
  
  # Fit gaussian process regression models for each response variable ----
  
  ## Myelin
  right_side <- '~ gp(age, k=10, c=5/4) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor'
  left_side <- paste0('myelin_glasser_v', region_id)
  
  ## Thickness
  # right_side <- '~ gp(age, k=10, c=5/4) + sex + site'
  # left_side <- paste0('thickness_glasser_v', region_id)
  
  mf_form <- as.formula(paste0(left_side, gsub('gp\\((\\w+).*\\)', '\\1', right_side)))
  form <- as.formula(paste0(left_side, right_side))
  
  d <- model.frame(formula = mf_form, data = hcpd_data)
  
  ## Create row index for merging subject IDs
  d_row_index <- as.numeric(row.names(d))
  d$subject_id <- hcpd_data$subject_id[d_row_index]
  agerange <- range(d$age)
  
  ## Identify column with response variable in model frame
  col_idx <- grep(paste0('myelin_glasser_v', region_id), colnames(d))
  
  ## Scale and center response variable
  d[, col_idx] <- scale(d[, col_idx], center = TRUE, scale = TRUE)
  
  ## Specify model output name
  model_name <- paste0('gp_fit_', left_side)
  
  # start_time <- Sys.time()
  
  ## Fit gaussian process regresion model with brms
  gp_fit <- brms::brm(brms::bf(form),
                      data=d,
                      chains = 4, cores = 4,
                      iter = 4500, warmup = 2000,
                      control = list(adapt_delta = .9999, max_treedepth = 20),
                      file = model_name, silent = FALSE)
  
  ## This is the function to use for prediction error (both deviation of the
  ## prediction from observed, and the error in that). That is, the output of this
  ## function gives us samples from y_obs - y_hat. The standard deviation of these
  ## samples (governed entirely by the sampling distribution of the expected mean,
  ## y_hat) is the standard error of the prediction (i.e, variance in the expected
  ## response), which, when squared, is what Marquand et al (2016) report as
  ## \sigma^2_{ij}. The normative differerence method also takes into account the
  ## error variance from the model's residuals, which is what Marquand et al
  ## (2016) refer to as "the variance learned from the normative distribution" and
  ## \sigma^2_{nj}. We can get the mean of the posterior distribution of
  ## \sigma^2_{nj} using the posterior_summary function.
  ##
  
  # 10-fold cross validation ----
  ## We can also use
  #?brms::kfold
  #?brms::kfold_predict()
  #?brms::predict.brmsfit
  #?brms::posterior_predict.brmsfit()
  #?brms::prepare_predictions
  
  #!!! For this next part to work properly, you need to make sure that the batch job gives you same number of CPUs as specified number of workers.
  
  #CHANGE THE NUMBER OF WORKERS TO MATCH YOUR MACHINE!
  if(!file.exists(paste0('kfold_gp_fit_', left_side, '.rds'))){
    future::plan(future::multiprocess, workers = 10)
    kfold1 <- brms::kfold(gp_fit, K = 10, save_fits = TRUE, chains = 4, nug = 1e-07)
    saveRDS(kfold1, paste0('kfold_gp_fit_', left_side, '.rds'))
  } else {
    kfold1 <- readRDS(paste0('kfold_gp_fit_', left_side, '.rds'))
  }
  
  # end_time <- Sys.time()
  # end_time - start_time

# }