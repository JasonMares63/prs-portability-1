#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=prep-PRS
#SBATCH -c 10
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=4gb

set -e

module load anaconda/3-4.2.0
source activate prs

# Basophil BMI (already done)
phenotypes=(DBP Eosinophil Hb Height Ht Lymphocyte MCH MCHC MCV Monocyte Neutrophil Platelet RBC SBP WBC)

# Threshold   0    1    2    3    4
thresholds=(5e-8 1e-6 1e-4 1e-3 1e-2)

# Clump GWAS results
for chromosome in $(seq 1 22);
do
  plink2 \
    --memory 35000 \
    --bgen /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/_001_ukb_imp_chr${chromosome}_v2.bgen ref-first \
    --sample /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/ukb_imp.sample \
    --keep data/ukb_populations/EUR_all.txt \
    --extract data/bbj/snps.txt \
    --make-just-pvar \
    --out data/ukb_bbj_pgen/${chromosome}

  python scripts/06a_find_duplicates.py \
    data/ukb_bbj_pgen/${chromosome}.pvar \
    --output data/ukb_bbj_pgen/${chromosome}_duplicates.txt

  rm data/ukb_bbj_pgen/${chromosome}.pvar

  plink2 \
    --memory 35000 \
    --bgen /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/_001_ukb_imp_chr${chromosome}_v2.bgen ref-first \
    --sample /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/ukb_imp.sample \
    --exclude data/ukb_bbj_pgen/${chromosome}_duplicates.txt \
    --make-bed \
    --out data/ukb_bbj_pgen/${chromosome}
done

for phenotype in "${phenotypes[@]}"
do
  # Clump GWAS results
  for chromosome in $(seq 1 22);
  do
    # Convert the GWAS output file from the Plink 2 format to Plink 1 `.assoc` format (removing duplicates)
    python scripts/06b_convert_plink2_glm_to_plink1.py \
      data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.linear \
      --output data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc

    # Clump GWAS results using LD from 1000 Genomes EUR super-population
    plink \
      --memory 35000 \
      --bfile data/ukb_bbj_pgen/${chromosome} \
      --clump data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc \
      --clump-p1 0.01 \
      --clump-p2 1 \
      --clump-r2 0.5 \
      --clump-kb 250

    # `--clump` automatically outputs plink.* files. Move and rename these
    mv plink.clumped data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.clumped
    mv plink.log data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.clumped.log
    mv plink.nosex data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.clumped.nosex
  done

  # Combine clumped SNPs across chromosomes
  head -n 1 data/gwas_results/${phenotype}.chr1.${phenotype}.clumped > \
    data/gwas_results/${phenotype}_combined.clumped
  tail -n+2 -q data/gwas_results/${phenotype}.chr*.${phenotype}.clumped >> \
    data/gwas_results/${phenotype}_combined.clumped

  # Create files of SNPs meeting several p-value thresholds. Files numbered 0-4, corresponding values above.
  for threshold in 0 1 2 3 4
  do
    # Further filter clumped SNPs using p-value thresholds (removes multiallelic SNPs)
    python scripts/06c_filter_snps_for_prs.py \
      data/gwas_results/${phenotype}_combined.clumped \
      --threshold ${thresholds[$threshold]} \
      --output data/prs/${phenotype}_threshold_${threshold}.txt
  done

  # Create combined GWAS result files for each phenotype
  python scripts/06d_combine_glm_threshold_4.py \
    /rigel/mfplab/users/mnz2108/prs-portability/data/gwas_results/${phenotype}.chr*.${phenotype}.glm.linear \
    --keep /rigel/mfplab/users/mnz2108/prs-portability/data/prs/${phenotype}_threshold_4.txt \
    --output /rigel/mfplab/users/mnz2108/prs-portability/data/gwas_results/${phenotype}_combined.glm.linear
done

######################################### Combine files to make a single PGEN ##########################################

# Combine all SNPs meeting the least extreme threshold for any phenotype
cat data/prs/*_threshold_4.txt | uniq > data/prs/all_prs_snps.txt

# Extract any SNPs meeting the least extreme criteria from all chromosomes
rm -f data/prs/prs_merged_list.txt
for chromosome in $(seq 1 22);
do
  plink \
    --bfile data/ukb_bbj_pgen/${chromosome} \
    --extract data/prs/all_prs_snps.txt \
    --exclude data/ukb_bbj_pgen/${chromosome}_duplicates.txt \
    --memory 35000 \
    --make-bed \
    --out data/prs/chr${chromosome}_temp

  printf "data/prs/chr%s_temp\n" $chromosome >> data/prs/prs_merged_list.txt
done

# Combine extracted SNPs across chromosomes into a single Plink 1 file of all
#  SNPs that meet the least extreme level of significance for any of the traits.
plink \
  --merge-list data/prs/prs_merged_list.txt \
  --make-bed \
  --memory 35000 \
  --out data/prs/merged

rm data/prs/chr*_temp.*
