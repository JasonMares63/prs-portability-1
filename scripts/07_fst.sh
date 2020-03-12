#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=fst
#SBATCH -c 2
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=6gb

set -e

module load anaconda/3-4.2.0
source activate prs

python scripts/07a_make_fst_clusters.py

plink \
  --bfile data/ukb_merged/merged \
  --within data/ukb_populations/fst_clusters.txt \
  --fst

mv plink.fst data/ukb_populations/
