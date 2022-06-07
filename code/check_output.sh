#!/bin/sh

gp_output_check=/ncf/hcp/data/analyses/myelin/Adversity_Project/results/output_check/missing_regional_gp_model_output.txt
rm $gp_output_check
touch $gp_output_check

kfold_output_check=/ncf/hcp/data/analyses/myelin/Adversity_Project/results/output_check/missing_regional_kfold_output.txt
rm $kfold_output_check
touch $kfold_output_check



for i in {1..360}
do
	## Check for normative modeling output
	outpath=/ncf/hcp/data/analyses/myelin/Adversity_Project/results/normative_modeling/regional/kfold_gp_fit_myelin_glasser_v${i}.rds

	if test -f "$outpath"; then
		echo "Good"
	else 
		echo "$i" >> $kfold_output_check
	fi
done

## How many missing models?
wc -l ${kfold_output_check}
