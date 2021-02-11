#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=merged_UKBB_filtered
#SBATCH -c 5
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=10gb

# Fail if any command fails
set -e

plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'

$plink \
  --merge-list data/ukb_merged/ukb_merged_list.txt \
  --make-bed \
  --out data/ukb_merged/merged

$plink2 \
  --bfile data/ukb_merged/merged \
  --make-pgen \
  --out data/ukb_merged/merged
