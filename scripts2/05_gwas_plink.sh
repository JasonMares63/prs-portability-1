#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=GWAS
#SBATCH -c 10
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=8gb

set -e

#Central directory version of plink not up to-date
plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'

for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC
do
  # Can't use Plink optimization for all phenotypes simultaneously because
  # each phenotype uses different samples.
  # The below hashtags refer to past GWAS, using genotype and UKBB
  for chromosome in $(seq 1 22);
  do
    $plink2\
      --bfile data/ukb_filtered/chr${chromosome}
#      --bgen /rigel/mfplab/projects/ukb_hakhamanesh/imputed/bgen_files/_001_ukb_imp_chr${chromosome}_v2.bgen ref-first \
#      --sample /rigel/mfplab/projects/ukb_hakhamanesh/imputed/bgen_files/ukb_imp.sample \
#      --extract /rigel/mfplab/projects/prs-portability/data/bbj/snps.txt \
#      --keep data/ukb_populations/${phenotype}.txt
      --keep data/ukb_populations/${phenotype}1.txt \
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
