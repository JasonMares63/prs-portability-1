#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=PCA
#SBATCH -c 3
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=8gb

# Fail if any command fails
set -e

#Central directory of plink software is not up to-date
plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454//plink/plink2'

# Compute PCA on the 1000 Genomes project filtered and merged genotypes
$plink2 \
  --bfile data/kgp_merged/merged \
  --freq counts \
  --pca allele-wts 20 \
  --out data/kgp_merged/merged

#  Apply loadings from PCA on 1000 Genomes to 1000 Genomes (for consistency)
$plink2 \
  --bfile data/kgp_merged/merged \
  --read-freq data/kgp_merged/merged.acount \
  --score data/kgp_merged/merged.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-25 \
  --out data/kgp_merged/projection

# Apply loadings from PCA computed on 1000 Genomes to the merged UK Biobank file
$plink2 \
  --bfile data/ukb_merged/merged \
  --read-freq data/kgp_merged/merged.acount \
  --score data/kgp_merged/merged.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-25 \
  --out data/ukb_merged/projection
