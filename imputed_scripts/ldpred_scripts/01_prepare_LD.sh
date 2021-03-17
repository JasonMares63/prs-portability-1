#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=ldpred_prepare
#SBATCH -c 2
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8gb

set -e

module load anaconda
module load R
source activate prs1

# Adapt bfile into LDpred objects
#Rscript ldpred_scripts/01a_create_ldpred_objects.R

# Format GWAS summary statistics (flipping betas to match A0/A1 of map)
#Rscript ldpred_scripts/01b_ldpred_gwas_sumstats.R

# Create LD reference panel correlation matrices
Rscript ldpred_scripts/01c_prepare_LD_matrices.R
