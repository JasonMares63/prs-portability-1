#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=Monocyte
#SBATCH -c 1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=10gb
set -e

module load anaconda
module load R
source activate prs1

Rscript ldpred_scripts/temp_scripts/Monocyte.R
rm ldpred_scripts/temp_scripts/Monocyte.R
