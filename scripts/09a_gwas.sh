#!/bin/bash
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=GWAS          # The job name.
#SBATCH -c 20                    # The number of cpu cores to use.
#SBATCH --time=1:00:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=5gb        # The memory the job will use per cpu core.

set -e

# module load anaconda/3-4.2.0
# source activate ../envs
#
# python scripts/09b_remove_psam_sex_column.py
# source deactivate

for phenotype in BMI DBP SBP Height
do
  plink2 \
    --pfile data/ukb_populations/EUR_all \
    --psam data/ukb_populations/EUR_all_nosex.psam \
    --keep data/ukb_populations/${phenotype}.txt \
    --pheno data/phenotypes/full_phenotypes.pheno \
    --pheno-name $phenotype \
    --covar data/ukb_filtered/covar_all_samples.covar \
    --covar-variance-standardize age sex_covar age_sq age_sex age_sq_sex PC1_AVG PC2_AVG PC3_AVG PC4_AVG PC5_AVG PC6_AVG PC7_AVG PC8_AVG PC9_AVG PC10_AVG PC11_AVG PC12_AVG PC13_AVG PC14_AVG PC15_AVG PC16_AVG PC17_AVG PC18_AVG PC19_AVG PC20_AVG \
    --require-pheno $phenotype \
    --require-covar age sex_covar age_sq age_sex age_sq_sex PC1_AVG PC2_AVG PC3_AVG PC4_AVG PC5_AVG PC6_AVG PC7_AVG PC8_AVG PC9_AVG PC10_AVG PC11_AVG PC12_AVG PC13_AVG PC14_AVG PC15_AVG PC16_AVG PC17_AVG PC18_AVG PC19_AVG PC20_AVG \
    --glm hide-covar \
    --out data/gwas_results/interaction
done
