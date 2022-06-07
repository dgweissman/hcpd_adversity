#!/bin/bash
#SBATCH -n 1
#SBATCH -p fasse
#SBATCH --signal=USR2
#SBATCH --time=4-00:00:00
#SBATCH --account=somerville_lab
#SBATCH --job-name=hcp_adversity_normative_modeling
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
# Outputs ----------------------------------
#SBATCH -o /ncf/hcp/data/analyses/myelin/containers/logs/%A_%a-%x.out

container=/n/home_fasse/gbaum/.fasseood/dev/ood-rstudio-versecmdstan/verse-cmdstan-ggseg-libs.simg

echo "singularity version $(singularity version)"

OVERLAY="$SCRATCH/LABS/somerville_lab/Users/gbaum/$(uuidgen).img"
srun -c 1 singularity overlay create --size 2512 "${OVERLAY}"

srun -c $SLURM_CPUS_PER_TASK singularity exec \
        --overlay ${OVERLAY} \
        ${container} \
        Rscript --no-save --no-restore /ncf/hcp/data/analyses/myelin/Adversity_Project/code/run_regional_normative_models.R
