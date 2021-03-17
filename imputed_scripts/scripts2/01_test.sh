#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=merge
#SBATCH -c 5
#SBATCH --time=15:00:00
#SBATCH --mem-per-cpu=16gb

set -e
plink="/rigel/mfplab/users/jm4454/plink/plink"
plink2="/rigel/mfplab/users/jm4454/plink/plink2"

$plink \
  --merge-list data/ukb_merged/ukb_merged_list.txt \
  --make-bed \
  --out data/ukb_merged/merged
