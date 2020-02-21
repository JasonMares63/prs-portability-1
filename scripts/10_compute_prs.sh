#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=PRS
#SBATCH -c 4
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=3gb

set -e

module load anaconda/3-4.2.0
source activate ../envs

for phenotype in BMI DBP SBP Height
do
  # Convert the GWAS output file from the Plink 2 format to Plink 1 `.assoc` format
  python scripts/10a_convert_plink2_glm_to_plink1.py \
    data/gwas_results/no_interaction.${phenotype}.glm.linear \
    --output data/gwas_results/no_interaction.${phenotype}.glm.assoc

  # Clump GWAS results using LD thresholds
  /rigel/mfplab/users/mnz2108/plink/plink \
    --bfile data/kgp_merged/merged \
    --clump data/gwas_results/no_interaction.${phenotype}.glm.assoc \
    --clump-p1 0.01 \
    --clump-p2 1 \
    --clump-r2 0.5 \
    --clump-kb 250

  # Clumping automatically outputs plink.* files. Move and rename these
  mv plink.clumped data/gwas_results/no_interaction.${phenotype}.clumped
  mv plink.log data/gwas_results/no_interaction.${phenotype}.log
  mv plink.nosex data/gwas_results/no_interaction.${phenotype}.nosex
done

# Threshold   0    1    2    3    4
thresholds=(5e-8 1e-6 1e-4 1e-3 1e-2)

for phenotype in BMI DBP SBP Height
do
  for threshold in 0 1 2 3 4
  do
    # Further filter clumped SNPs using p-value thresholds
    python scripts/10b_filter_snps_for_prs.py \
      data/gwas_results/no_interaction.${phenotype}.clumped \
      --threshold ${thresholds[$threshold]} \
      --output data/gwas_results/PRS_${phenotype}_threshold_${threshold}.txt

    # Score each individual using the SNPs below a p-value threshold and GWAS betas
    /rigel/mfplab/users/mnz2108/plink/plink2 \
      --pfile data/ukb_filtered/merged \
      --extract data/gwas_results/PRS_${phenotype}_threshold_${threshold}.txt \
      --score data/gwas_results/no_interaction.${phenotype}.glm.linear 3 6 9 header no-mean-imputation \
      --out data/gwas_results/PRS_${phenotype}_${threshold}_scores
  done
done
