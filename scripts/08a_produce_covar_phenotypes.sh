#!/bin/bash
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=covarPhenotype      # The job name.
#SBATCH -c 2                     # The number of cpu cores to use.
#SBATCH --time=60:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=6gb        # The memory the job will use per cpu core.

set -e

module load anaconda/3-4.2.0
source activate geno

# Creates covariate file with age, sex, age*sex, age^2, age^2 * sex and PC1, ..., PC20
python scripts/08b_create_covariates.py

# Creates data/phenotypes/full_phenotypes.pheno, combining all phenotypes for GWAS 
Rscript scripts/08c_create_phenotypes_file.R

