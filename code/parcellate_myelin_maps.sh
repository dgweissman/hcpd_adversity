#!/bin/sh

## Load Connectome Workbench module
module load connectome-workbench/1.3.2-fasrc01

## Define input data
cifti_in=/ncf/hcp/data/analyses/myelin/Adversity_Project/data/myelin_maps/HCD1320_Winter2022/HCD1320_Winter2022.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR.dscalar.nii

# cifti_in=/ncf/hcp/data/analyses/myelin/parcellations/SensorimotorAssociation.Axis.Glasser360.pscalar.nii

# Parcellations
glasser_atlas=/ncf/hcp/data/analyses/myelin/parcellations/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii
yeo7_atlas=/ncf/hcp/data/analyses/myelin/parcellations/RSN-networks.32k_fs_LR.dlabel.nii
cole_atlas=/ncf/hcp/data/analyses/myelin/parcellations/CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_netassignments_LR.dlabel.nii
schaefer_atlas=/ncf/hcp/data/analyses/myelin/parcellations/Schaefer2018_400Parcels_7Networks_order.dlabel.nii

# Define Output
outdir=/ncf/hcp/data/analyses/myelin/Adversity_Project/data/myelin_maps/

cifti_out=${outdir}/HCD1320_Winter2022.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Cole-Anticevic-Networks.pscalar.nii

text_out=${outdir}/HCD1320_Winter2022.All.MyelinMap_IndPseudoCorr_MSMAll.32k_fs_LR_Cole-Anticevic-Networks.txt

## Extract mean myelin values for each parcel in brain atlas
wb_command -cifti-parcellate ${cifti_in} ${cole_atlas} COLUMN ${cifti_out} -method MEAN

## Output parcel values as text file
wb_command -cifti-convert -to-text ${cifti_out} ${text_out} 

