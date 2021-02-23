
import os

variable = """#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name={}
#SBATCH -c 1
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=10gb
set -e
plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'
phenotype={}
"""

fixed = """# Create PRS score files for several p-value thresholds. Files numbered 0-4, corresponding values above.
#for i in 6 $(seq 18 28);
for i in 6 28;
do
$plink2 \
--bfile data/ukb_merged/merged \
--score data/LDpred/prs/${phenotype}_LDpred_sumstats.tsv 2 5 ${i} header no-mean-imputation \
--memory 10000 \
--out data/LDpred/prs/${phenotype}_LDpred_scores_${i}
done

"""


def main(phenotype):
    return variable.format(phenotype, phenotype) + fixed


if __name__ == "__main__":
    phenotypes = ['BMI','Height','RBC', 'Platelet', 'MCV', 'Monocyte', 'WBC', 'MCH', 'Eosinophil', 'Lymphocyte']
    for phenotype in phenotypes:
        with open(f'ldpred_scripts/temp_scripts/{phenotype}.sh', 'w') as f:
            f.write(main(phenotype))
        os.system(f'sbatch ldpred_scripts/temp_scripts/{phenotype}.sh')
