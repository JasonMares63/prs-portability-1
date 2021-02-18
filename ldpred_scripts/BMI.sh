#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=BMI
#SBATCH -c 1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=10gb
set -e

module load anaconda
module load R
source activate prs1

Rscript ldpred_scripts/temp_scripts/BMI.R
rm ldpred_scripts/temp_scripts/BMI.R
