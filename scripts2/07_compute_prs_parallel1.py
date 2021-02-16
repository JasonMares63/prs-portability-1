import os


variable = """#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name={}
#SBATCH -c 1
#SBATCH --time=30:00
#SBATCH --mem-per-cpu=6gb

set -e

plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'

phenotype={}

"""

fixed = """# Create PRS score files for several p-value thresholds. Files numbered 0-4, corresponding values above.
for threshold in 0 1 2 3 4
do
  # Score each individual using the SNPs below a p-value threshold and GWAS betas
  $plink2 --bfile data/prs/test_UKBB_imp_BBJ_GWAS_1KG_EUR_LD_merged --extract data/prs/${phenotype}_threshold_${threshold}.txt \
    --score data/gwas_results/${phenotype}_combined.glm.linear 3 6 9 header no-mean-imputation \
    --memory 10000 \
    --out data/prs/${phenotype}_${threshold}_scores
done
"""


def main(phenotype):
    return variable.format(phenotype, phenotype) + fixed


if __name__ == "__main__":
    phenotypes = ['BMI','Basophil','DBP', 'Eosinophil', 'Hb', 'Height', 'Ht', 'Lymphocyte', 'MCH', 'MCHC', 'MCV', 'Monocyte', 'Neutrophil', 'Platelet', 'RBC', 'SBP', 'WBC']
    phenotypes1 = ['Height', 'BMI' ,'RBC', 'Platelet', 'MCV', 'Monocyte', 'WBC', 'MCH', 'Eosinophil', 'Lymphocyte']
    for phenotype in phenotypes1:
        with open(f'{phenotype}.sh', 'w') as f:
            f.write(main(phenotype))
        os.system(f'sbatch {phenotype}.sh')
