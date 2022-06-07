rm(list = ls())

##############################
# Load relevant Libraries ----
##############################
require(ggplot2)
require(stringr)

####################################################
## Read in demographics and data quality measures ##
####################################################
# harms_QA <- read.csv("/ncf/hcp/data/analyses/myelin/data/subject_demographics/freesurfer_qc_20191210.csv", header=TRUE)
# hcpd_qa <- subset(harms_QA, harms_QA$Project =="HCD")

hcpd_data <- read.csv("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/subject_demographics/hcpd_n1320_demographics_wlogINR.csv")
colnames(hcpd_data)[6] 
colnames(hcpd_data)[6] <- "sex"
## Set variables
hcpd_data$sex<- as.factor(hcpd_data$sex)
hcpd_data$site <- as.factor(hcpd_data$site)
# hcpd_data$Scanner  <- as.factor(hcpd_data$Scanner)
hcpd_data$SES_RLVL <- as.factor(hcpd_data$SES_RLVL)

# Create age rounded down to integer as factor
hcpd_data$age_int <- as.factor(round(hcpd_data$age))
hcpd_data$age_floor <- as.factor(floor(hcpd_data$age))

## Subset by first timepoint and 8-21 years old
subject_index <- which(hcpd_data$timepoint=="V1" & hcpd_data$age >= 8)
length(subject_index)

# hcpd_data <- hcpd_data[subject_index, ]

hcpd_data <- subset(hcpd_data, hcpd_data$timepoint=="V1" & hcpd_data$age >= 8)
dim(hcpd_data)


## Look at age distribution
hist(hcpd_data$age, col="slategray3")
table(hcpd_data$sex)

# hcpd_data$Subj.ID <- hcpd_data$Subject

######################################
## READ IN GLASSER'S CORRECTED MAPS ##
######################################
library(cifti)
library(data.table)
HCD_myelin_maps_cifti <- read_cifti("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/myelin_maps/HCD1320_Winter2022/HCD1320_Winter2022.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR.dscalar.nii")

myelin_maps <- HCD_myelin_maps_cifti$data
myelin_maps <- transpose(as.data.frame(myelin_maps))
mean_myelin_map <- unlist(lapply(myelin_maps, mean))

map_names <- HCD_myelin_maps_cifti$NamedMap
id <- substr(map_names$map_names, 1, 10)
subject_id <- substr(map_names$map_names, 1, 16)

#################################################
## Read in Glasser-MMP Parcellated Myelin Maps ##
#################################################
GlasserMMP_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/myelin_maps/HCD1320_Winter2022.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Glasser-MMP.txt")
GlasserMMP_transmitCorr_myelinMaps <- as.data.frame(t(as.matrix(GlasserMMP_transmitCorr_myelinMaps)))

## Re-name columns for regional myelin
for(i in 1:360) {
  colnames(GlasserMMP_transmitCorr_myelinMaps)[i] <- paste("myelin_glasser_v", i, sep = "")
}

GlasserMMP_transmitCorr_myelinMaps$subject_id <- subject_id
GlasserMMP_transmitCorr_myelinMaps$subject_id <- noquote(GlasserMMP_transmitCorr_myelinMaps$subject_id)


## Retain subjects from timepoint 1 and ages 8-21
# GlasserMMP_transmitCorr_myelinMaps <- GlasserMMP_transmitCorr_myelinMaps[subject_index,]
# hcpd_data <- cbind(hcpd_data, GlasserMMP_transmitCorr_myelinMaps)

hcpd_data <- merge( hcpd_data, GlasserMMP_transmitCorr_myelinMaps,by="subject_id")
dim(hcpd_data)


########################################
## Read in Cole-Anticevic Myelin Maps ##
########################################
ColeAnt_12network_transmitCorr_myelinMaps <- read.table("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/myelin_maps/HCD1320_Winter2022.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Cole-Anticevic-Networks.txt")

ColeAnt_12network_transmitCorr_myelinMaps <- as.data.frame(t(as.matrix(ColeAnt_12network_transmitCorr_myelinMaps)))
dim(ColeAnt_12network_transmitCorr_myelinMaps)

## Re-name columns for Yeo7 network myelin
colnames(ColeAnt_12network_transmitCorr_myelinMaps) <- c("Visual1_myelin", "Visual2_myelin", "Somatomotor_myelin", "Cingulo_Opercular_myelin", "Dorsal_attention_myelin", "Language_myelin", "Frontoparietal_myelin", "Auditory_myelin", "Default_myelin", "Posterior_Multimodal_myelin", "Ventral_Multimodal_myelin", "Orbito_Affective_myelin")

ColeAnt_12network_transmitCorr_myelinMaps$subject_id <- subject_id

hcpd_data <- merge(hcpd_data, ColeAnt_12network_transmitCorr_myelinMaps, by="subject_id")
dim(hcpd_data)

##################################
## Read in Nuissance Covariates ##
##################################

glasser_covars <- read.csv("/ncf/hcp/data/analyses/myelin/Adversity_Project/data/myelin_maps/HCD1320_Winter2022/n1320_Covariates.csv")
glasser_covars$subject_id <- subject_id

## Merge covariates with df
hcpd_data <- merge(hcpd_data, glasser_covars, by="subject_id")

names(hcpd_data)
dim(hcpd_data)

############################
## Mean Wholebrain Myelin ##
############################ 
nreg <- 360
myelin_col_index <- grep("myelin_glasser_v", colnames(hcpd_data))

# across regions (within subjects)
tmp_myelin <-as.data.frame(hcpd_data[myelin_col_index])
mean_wholebrain_T1wT2w <- rowMeans(tmp_myelin)
hcpd_data$mean_wholebrain_T1wT2w <- mean_wholebrain_T1wT2w
hist(hcpd_data$mean_wholebrain_T1wT2w, col="slategray3")

## Create Threshold Based on Wholebrain T1wT2w ##
high_wholebrain_thresh <- mean(hcpd_data$mean_wholebrain_T1wT2w) + (4*sd(hcpd_data$mean_wholebrain_T1wT2w))
low_wholebrain_thresh <- mean(hcpd_data$mean_wholebrain_T1wT2w) - (4*sd(hcpd_data$mean_wholebrain_T1wT2w))

## Look at subjects with outlier data
myelin_thresh_idx <- which(hcpd_data$mean_wholebrain_T1wT2w > high_wholebrain_thresh | hcpd_data$mean_wholebrain_T1wT2w < low_wholebrain_thresh)
length(myelin_thresh_idx)
hcpd_data$subject_id[myelin_thresh_idx]

#######################################################
## REMOVE SUBEJECTS WITH WHOLEBRAIN T1w/T2w OUTLIERS ##
#######################################################
# hcpd_data <- hcpd_data[-myelin_thresh_idx,]
# dim(hcpd_data)

####################################
## Read in cortical thickness maps ##
#####################################
nsub <- dim(hcpd_data)[1]
nreg <- 360
glasser_thickness_maps <- as.data.frame(array(rep(NA, nsub*nreg), dim=c(nsub, nreg)))

# for(i in 1:nsub) {
#   sub_id <- hcpd_data$subject_id[i]
#   glasser_thickness_path <- paste0("/users/gbaum/myelin_project/intradb_output/n632_thickness/", sub_id, "_thickness.32k_fs_LR_Glasser-MMP.txt")
#   tmp_thickness <-  read.table(glasser_thickness_path, header=FALSE) 
#   glasser_thickness_maps[i,] <- t(tmp_thickness)
# }

## Re-name columns for regional myelin
# for(i in 1:360) {
#   colnames(glasser_thickness_maps)[i] <- paste("thickness_glasser_v", i, sep = "")
# }

## Merge with hcpd_data
# hcpd_data <- cbind(hcpd_data, glasser_thickness_maps)

#######################################
# Visualize Sample Characteristics ----
#######################################

## Graphics setup
apal <- paste0('#', c('000000', 'EAE3D8', 'FFFFFF', 'FFB838', '486D87'))

gbtheme <- theme_classic() +  
  theme(text = element_text(size = 36),
        panel.background = element_rect(fill = apal[[3]], size = 0, color = apal[[2]]),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = apal[[2]], size = 0),
        strip.text = element_text(color = '#222222'),
        axis.text =  element_text(color = apal[[1]]), axis.title = element_text(color = apal[[1]]),
        axis.ticks.length=unit(0.5, "cm"))

## Age by Sex: Proportion
age_sex_plot <- ggplot(hcpd_data, aes(x=age_floor, fill=sex, col=sex)) + geom_bar(position = "fill", alpha=0.7) + scale_y_continuous(labels = scales::percent)
age_sex_plot + gbtheme + scale_fill_manual(values = c("#424C6D",  "#CE7B5B")) + scale_colour_manual(values = c("#424C6D",  "#CE7B5B"))  

## Age by Sex: Count
age_sex_plot <- ggplot(hcpd_data, aes(x=age_floor, fill=sex, col=sex)) + geom_bar(alpha=0.7)
age_sex_plot + gbtheme + scale_fill_manual(values = c("#424C6D",  "#CE7B5B")) + scale_colour_manual(values = c("#424C6D",  "#CE7B5B"))  

## Age by SES: Proportion
age_ses_plot <- ggplot(hcpd_data, aes(x=age_floor, fill=SES_RLVL, col=SES_RLVL)) + geom_bar(position = "fill", alpha=0.7) + scale_y_continuous(labels = scales::percent)
age_ses_plot + gbtheme + scale_fill_manual(values = c("#424C6D", "#CE7B5B", "#FFD378", "#FFFFFF")) + scale_colour_manual(values = c("#424C6D", "#CE7B5B", "#FFD378", "#FFFFFF"))  

## Age by SES: Count
age_ses_plot <- ggplot(hcpd_data, aes(x=age_floor, fill=SES_RLVL, col=SES_RLVL)) + geom_bar(alpha=0.7)
age_ses_plot + gbtheme + scale_fill_manual(values = c("#424C6D", "#CE7B5B", "#FFD378", "#FFFFFF")) + scale_colour_manual(values = c("#424C6D", "#CE7B5B", "#FFD378", "#FFFFFF"))  

library(mgcv)
library(visreg)
## Test GAM using mean cortical myelin
test_gam <- gam(mean_wholebrain_T1wT2w ~  logINR + s(age) + sex + site + Reference.Voltage + Mean.Pseudo.Transmit.Map + T2..Dropout.Threshold + Smoothing.FWHM.mm + Pseudotransmit.Reference.Value + Correction.Equation.Slope + Corrected.CSF.Regressor,data=hcpd_data)
visreg(test_gam, "logINR")
visreg(test_gam, "age")


## Export data as csv
write.csv(hcpd_data, "/ncf/hcp/data/analyses/myelin/Adversity_Project/data/subject_demographics/n925_hcpd_myelin_May2022.csv", row.names=FALSE)
write.table(hcpd_data$subject_id, "/ncf/hcp/data/analyses/myelin/Adversity_Project/data/subject_demographics/n925_hcpd_myelin_subject_list.txt", row.names=FALSE, col.names=FALSE, quote = FALSE)
