#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=ldpred_adj
#SBATCH -c 1
#SBATCH --time=01-00:00:00
#SBATCH --mem-per-cpu=6gb

set -e

module load anaconda
module load R
source activate prs1

Rscript ldpred_scripts/03a_adjust_ldpred_auto_results.R
python ldpred_scripts/03b_ldpred_scores_batch.py

