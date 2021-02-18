#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=BMI
#SBATCH -c 1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=10gb


set -e
plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'

phenotype=BMI


#$plink2 \
#--bfile data/ukb_merged/merged \
#--keep data/ukb_populations/AFR_all.txt data/ukb_populations/AMR_all.txt data/ukb_populations/EAS_all.txt data/ukb_populations/EUR_test.txt data/ukb_populations/SAS_all.txt \
#--memory 30000 \
#--make-bed \
#--out data/ukb_merged/non_EUR_train

for i in 6 $(seq 24 28);
do
$plink2 \
--bfile data/ukb_merged/non_EUR_train \
--score data/LDpred2/prs/${phenotype}_LDpred_sumstats.tsv 2 5 ${i} header no-mean-imputation \
--memory 10000 \
--out data/LDpred2/prs/${phenotype}_LDpred_scores_${i}
done 

 

