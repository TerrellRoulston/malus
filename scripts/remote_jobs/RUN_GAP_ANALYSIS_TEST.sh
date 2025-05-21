#!/bin/bash
#SBATCH --job-name=malus_gap_array
#SBATCH --array=0           # <-- Only run index 0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --mail-user=terrell.roulston@acadiau.ca
#SBATCH --mail-type=ALL
#SBATCH --output=logs/malus_gap_%A_%a.out
#SBATCH --error=logs/malus_gap_%A_%a.err

# Load R module and libraries
module load StdEnv/2023 r/4.4.0
export R_LIBS=/project/6074193/mig_lab/bin/RPackages

# Set working directory to your project folder
cd /project/6074193/mig_lab/malus_gap/

# Run the R script and pass the SLURM array index
Rscript scripts/run_gap_analysis.R $SLURM_ARRAY_TASK_ID
