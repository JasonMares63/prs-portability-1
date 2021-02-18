#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=Platelet
#SBATCH -c 1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=10gb
set -e

module load anaconda
module load R
source activate prs1

Rscript ldpred_scripts/temp_scripts/Platelet.R
rm ldpred_scripts/temp_scripts/Platelet.R
