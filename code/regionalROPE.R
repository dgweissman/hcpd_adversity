library(brms)
region.logINR.ROPE<-data.frame(matrix(NA,nrow = 360,ncol = 11))
colnames(region.logINR.ROPE)[1]<-"region"
colnames(region.logINR.ROPE)[2:11]<-paste("pbeta",substr(seq(.01,.1,by=.01),2,4),sep="")
region.p_grade.ROPE<-region.logINR.ROPE
a<-commandArgs(trailingOnly=TRUE)
hcpd_data <- read.csv("n925_hcpd_myelin_AS.csv")
roi <- as.numeric(a)
error_index <- which(hcpd_data$p_grade_cat==77)
hcpd_data$subject_id[error_index]
hcpd_data$p_grade_cat[error_index] <- NA
hcpd_data$p_grade[error_index] <- NA
hcpd_data$p_grade_cat <- droplevels.factor(hcpd_data$p_grade_cat)

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
region<-colnames(hcpd_data[,myelin_col_index])[roi]

glasser_gam_results <- NULL
mod <- NULL

## regional logINR ROPE Test 
prior1 <- c(set_prior("normal(0, 10)", class = "b"), set_prior("normal(0, 10)", class = "Intercept"))
m <- NULL
full_covars=" ~  logINR + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
m <- as.formula(paste(region, full_covars, sep=""))
region.logINR.brm <- brm(formula=m,prior = prior1, data = hcpd_data, family = gaussian(link = "identity"), warmup = 1e3, iter = 1.5e4, thin = 5, chains = 4, cores = 4, seed = "123", control = list(adapt_delta = 0.999, max_treedepth = 20))
stdratio<-sd(hcpd_data$logINR,na.rm = T)/sd(hcpd_data[,region],na.rm=T)
region.logINR.ROPE<-read.csv("region.logINR.ROPE.csv")
region.logINR.ROPE[a,1]<-region
samples.logINR.brm <- posterior_samples(region.logINR.brm, "logINR")
for (b in 1:10){
  region.logINR.ROPE[a,b+1]<-mean(samples.logINR.brm$b_logINR*stdratio > -.01*b & samples.logINR.brm$b_logINR*stdratio < .01*b)
  }
#write.csv(region.logINR.ROPE,"region.logINR.ROPE.csv",row.names=F)
## regional p_grade ROPE Test 
m <- NULL
full_covars=" ~  p_grade + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor"   
m <- as.formula(paste(region, full_covars, sep=""))
region.p_grade.brm <- brm(formula=m,prior = prior1, data = hcpd_data, family = gaussian(link = "identity"), warmup = 1e3, iter = 1.5e4, thin = 5, chains = 4, cores = 4, seed = "123", control = list(adapt_delta = 0.999, max_treedepth = 20))
stdratio<-sd(hcpd_data$p_grade,na.rm = T)/sd(hcpd_data[,region],na.rm=T)
region.p_grade.ROPE<-read.csv("region.p_grade.ROPE.csv")
region.p_grade.ROPE[a,1]<-region
samples.p_grade.brm <- posterior_samples(region.p_grade.brm, "p_grade")
for (b in 1:10){
  region.p_grade.ROPE[a,b+1]<-mean(samples.p_grade.brm$b_p_grade*stdratio > -.01*b & samples.p_grade.brm$b_p_grade*stdratio < .01*b)
}
write.csv(region.p_grade.ROPE,"region.p_grade.ROPE.csv",row.names=F)