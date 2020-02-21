#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=covariates_phenotypes
#SBATCH -c 2
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=6gb

set -e

module load anaconda/3-4.2.0
source activate ../envs

# Creates covariate file with age, sex, age*sex, age^2, age^2 * sex and PC1, ..., PC20
python scripts/04a_create_covariates.py

# Creates data/phenotypes/full_phenotypes.pheno, combining all phenotypes for GWAS
# and inverse rank normal transforming phenotypes
Rscript scripts/04b_create_phenotypes_file.R
