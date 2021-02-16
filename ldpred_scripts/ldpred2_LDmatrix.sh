#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=ld_matrix
#SBATCH -c 6
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=12gb

module load anaconda
source activate prs1

Rscript LDmatrix.R

