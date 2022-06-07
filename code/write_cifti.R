########################################################
## Define CIFTI template made in Connectome Workbench ##
########################################################
glasser_template.path <- "/ncf/hcp/data/analyses/myelin/parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_template.pscalar.nii"
schaefer_template.path <- "/ncf/hcp/data/analyses/myelin/parcellations/S1200_Schaefer400x7_template.pscalar.nii"
# glasser_template.path <- "/ncf/hcp/data/analyses/myelin/parcellations/S1200_MMPtemplate.pscalar.nii"

#####################################################
##  Make a template pscalar.nii from the MMP atlas ##
#####################################################
system("module load connectome-workbench/1.3.2-fasrc01; wb_command -cifti-parcellate /ncf/hcp/data/analyses/myelin/parcellations/S1200.MyelinMap_BC_MSMAll.32k_fs_LR.dscalar.nii /ncf/hcp/data/analyses/myelin/parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii COLUMN /ncf/hcp/data/analyses/myelin/parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_template.pscalar.nii");

#############################################################
## Function to write CIFTI using a vector of parcel values ##
#############################################################
write_cifti <- function(template_path, brainVar_path, cifti_output_path) {
  
  system(paste0("module load connectome-workbench/1.3.2-fasrc01; ", "wb_command -cifti-convert -from-text ", brainVar_path, " ", template_path, " ", cifti_output_path))
  return(print(cifti_output_path))
}

#############################
## Edit input/output paths ##
#############################
brainVar.path <- "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/n925_logINR_normDiff_myelin_tstat.txt"
cifti.output.path <- "/ncf/hcp/data/analyses/myelin/Adversity_Project/results/cifti/n925_logINR_normDiff_myelin_tstat.pscalar.nii" 
glasser_template.path <- "/ncf/hcp/data/analyses/myelin/parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_template.pscalar.nii"

write_cifti(glasser_template.path, brainVar.path, cifti.output.path)
