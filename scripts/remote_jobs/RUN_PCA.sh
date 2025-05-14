#!/bin/bash
#SBATCH --job-name=malus_equiv_array
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G             
#SBATCH --time=02:00:00
#SBATCH --mail-user=terrell.roulston@acadiau.ca
#SBATCH --mail-type=ALL
#SBATCH --output=logs/malus_equiv_%A_%a.out
#SBATCH --error=logs/malus_equiv_%A_%a.err

# Load R Module and Packages
module load StdEnv/2023 r/4.4.0
export R_LIBS=/project/6074193/mig_lab/bin/RPackages

# Set dir to Malus PCA subfolder
cd /project/6074193/mig_lab/malus_pca/

# Run R Script
# Slurm Index [0,1,2] determines species pairs inside R script
Rscript scripts/run_malus_pca_equiv.R $SLURM_ARRAY_TASK_ID