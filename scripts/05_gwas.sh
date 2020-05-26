#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=GWAS
#SBATCH -c 20
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=2gb

set -e

for phenotype in Basophil BMI DBP Eosinophil Hb Height Ht Lymphocyte MCH MCHC MCV Monocyte Neutrophil Platelet RBC SBP WBC
do
  # Can't use Plink optimization for all phenotypes simultaneously because
  # each phenotype uses different samples.
  for chromosome in $(seq 1 22);
  do
    plink2 \
      --bgen /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/_001_ukb_imp_chr${chromosome}_v2.bgen ref-first \
      --sample /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/ukb_imp.sample \
      --extract data/bbj/snps.txt \
      --keep data/ukb_populations/${phenotype}.txt \
      --pheno data/phenotypes/full_phenotypes.pheno \
      --pheno-name $phenotype \
      --require-pheno $phenotype \
      --covar data/ukb_merged/covar_all_samples.covar \
      --covar-name age sex_covar age_sq age_sex age_sq_sex $(printf "PC%i_AVG " $(seq 1 20)) \
      --require-covar age sex_covar age_sq age_sex age_sq_sex $(printf "PC%i_AVG " $(seq 1 20)) \
      --vif 100000 \
      --memory 35000 \
      --glm hide-covar \
      --out data/gwas_results/${phenotype}.chr${chromosome}
  done
done
