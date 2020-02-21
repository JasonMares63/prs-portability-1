#!/bin/sh
#
#SBATCH --account=mfplab
#SBATCH --job-name=pca_apply
#SBATCH -c 10
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8gb


# Fail if any command fails
set -e

# Apply loadings from PCA computed on 1000 Genomes to the merged UK Biobank Plink 2 file
/rigel/mfplab/users/mnz2108/plink/plink2 \
  --pfile data/ukb_merged/merged \
  --read-freq data/kgp_merged/merged.acount \
  --score data/kgp_merged/merged.eigenvec.allele 2 5 header-read no-mean-imputation \
  --score-col-nums 6-25 \
  --out data/ukb_merged/projection

#  Apply loadings from PCA on 1000 Genomes to 1000 Genomes (for consistent method)
/rigel/mfplab/users/mnz2108/plink/plink2 \
  --pfile data/kgp_merged/merged \
  --read-freq data/kgp_merged/merged.acount \
  --score data/kgp_merged/merged.eigenvec.allele 2 5 header-read no-mean-imputation \
  --score-col-nums 6-25 \
  --out data/kgp_merged/projection
