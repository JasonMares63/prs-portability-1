import os


template = """#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=PRS-{}
#SBATCH -c 2
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=8gb

set -e

plink="/rigel/mfplab/users/jm4454/plink/plink"
plink2="/rigel/mfplab/users/jm4454/plink/plink2"

module load anaconda/3-4.2.0
source activate prs

# Threshold   0    1    2    3    4
thresholds=(5e-8 1e-6 1e-4 1e-3 1e-2)

phenotype={}
"""

remainder = """
for chromosome in $(seq 1 22);
do
  $plink2 \
    --bgen /rigel/mfplab/projects/ukb_hakhamanesh/imputed/bgen_files/_001_ukb_imp_chr${chromosome}_v2.bgen ref-first \
    --sample /rigel/mfplab/projects/ukb_hakhamanesh/imputed/bgen_files/ukb_imp.sample \
    --extract data/prs/${phenotype}_chr${chromosome}_threshold_4.txt \
    --memory 10000 \
    --make-pgen \
    --out data/prs/${phenotype}_chr${chromosome}_temp

  # Create PRS score files for several p-value thresholds. Files numbered 0-4, corresponding values above.
  for threshold in 0 1 2 3 4
  do
    # Further filter clumped SNPs using p-value thresholds
    python scripts/06b_filter_snps_for_prs.py \
      data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.clumped \
      --threshold ${thresholds[$threshold]} \
      --output data/prs/${phenotype}_chr${chromosome}_threshold_${threshold}.txt
  
    # Score each individual using the SNPs below a p-value threshold and GWAS betas
    $plink2 \
      --pfile data/prs/${phenotype}_chr${chromosome}_temp \
      --extract data/prs/${phenotype}_chr${chromosome}_threshold_${threshold}.txt \
      --score data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.linear 3 6 9 header no-mean-imputation ignore-dup-ids \
      --memory 10000 \
      --out data/prs/${phenotype}_chr${chromosome}_${threshold}_scores || true
  done
  rm data/prs/${phenotype}_chr${chromosome}_temp*
done
"""

phenotypes = ['BMI', 'Basophil', 'DBP', 'Eosinophil', 'Hb', 'Height', 'Ht', 'Lymphocyte', 
              'MCH', 'MCHC', 'MCV', 'Monocyte', 'Neutrophil', 'Platelet', 'RBC', 'SBP', 'WBC']

contents = {f'prs/{phenotype}.sh': template.format(phenotype, phenotype) + remainder 
            for phenotype in phenotypes}

os.system('mkdir prs/')

for filename, text in contents.items():
    with open(filename, 'w') as f:
        f.write(text)
    os.system(f'sbatch {filename}')

