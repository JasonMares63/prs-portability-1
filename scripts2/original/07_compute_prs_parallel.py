import os


variable = """#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name={}
#SBATCH -c 1
#SBATCH --time=30:00
#SBATCH --mem-per-cpu=12gb

set -e

plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'

phenotype={}

thresholds=(5e-8 1e-6 1e-4 1e-3 1e-2)

"""

fixed = """# Create PRS score files for several p-value thresholds. Files numbered 0-4, corresponding values above.
for threshold in 0 1 2 3 4
do
  # Score each individual using the SNPs below a p-value threshold and GWAS betas
  $plink \
    --bfile data/prs/merged \
    --extract data/prs/${phenotype}_threshold_${threshold}.txt \
    --score data/gwas_results/${phenotype}_combined.glm.linear 3 6 9 header no-mean-imputation \
    --memory 10000 \
    --out data/prs/${phenotype}_${threshold}_scores
done
"""


def main(phenotype):
    return variable.format(phenotype, phenotype) + fixed


if __name__ == "__main__":
    phenotypes = ['Basophil', 'BMI', 'DBP', 'Eosinophil', 'Hb', 'Height', 'Ht',
                  'Lymphocyte', 'MCH', 'MCHC', 'MCV', 'Monocyte', 'Neutrophil',
                  'Platelet', 'RBC', 'SBP', 'WBC']

    for phenotype in phenotypes:
        with open(f'prs/{phenotype}.sh', 'w') as f:
            f.write(main(phenotype))
        os.system(f'sbatch prs/{phenotype}.sh')
