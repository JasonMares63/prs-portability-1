#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=prep-PRS-ukb-bbj-gwas-ukb-ld
#SBATCH -c 10
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4gb

set -e

#Central directory of plink not up to-date
plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'
pipe='/rigel/mfplab/projects/prs-portability'
module load anaconda/3-4.2.0

source activate prs1

phenotypes=(Basophil BMI DBP Eosinophil Hb Height Ht Lymphocyte MCH MCHC MCV Monocyte Neutrophil Platelet RBC SBP WBC)
#phenotypes=(BMI Eosinophil Height Lymphocyte MCH MCV Monocyte Platelet RBC WBC)

# Threshold   0    1    2    3    4
thresholds=(5e-8 1e-6 1e-4 1e-3 1e-2)


###################### Already done
# Clump GWAS results
#for chromosome in $(seq 1 22);
#do
#  $plink2\
#    --memory 35000 \
#    --bfile data/kgp_filtered/chr${chromosome} \
#    --keep data/kgp_populations/EUR_all.txt \
#    --make-pgen \
#    --out data/kgp_populations/EUR_${chromosome}

#  python scripts2/06a_find_duplicates.py \
#    data/kgp_populations/EUR_${chromosome}.pvar \
#    --output data/ukb_bbj_pgen/${chromosome}_duplicates.txt

#  $plink2\
#    --memory 35000 \
#    --pfile data/kgp_populations/EUR_${chromosome} \
#    --exclude data/ukb_bbj_pgen/${chromosome}_duplicates.txt \
#    --make-bed \
#    --out data/kgp_populations/EUR_${chromosome}

#  rm data/kgp_populations/EUR_${chromosome}.p*

######################

for phenotype in "${phenotypes[@]}"
do
  # Clump GWAS results
  for chromosome in $(seq 1 22);
  do
    # Convert the GWAS output file from the Plink 2 format to Plink 1 `.assoc` format (removing duplicates)
#    python scripts/06b_convert_plink2_glm_to_plink1.py \
#      data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.linear \
#      --output data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc

    # Clump GWAS results using LD from 1000 Genomes EUR super-population
    $plink\
      --memory 35000 \
      --bfile data/kgp_populations/EUR_${chromosome} \
      --clump $pipe/data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc \
      --clump-p1 0.01 \
      --clump-p2 1 \
      --clump-r2 0.5 \
      --clump-kb 250 \
      --out data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}
  done

  # Combine clumped SNPs across chromosomes
  head -n 1 data/gwas_results/${phenotype}.chr1.${phenotype}.clumped > \
    data/gwas_results/${phenotype}_combined.clumped
  tail -n +2 -q data/gwas_results/${phenotype}.chr*.${phenotype}.clumped >> \
    data/gwas_results/${phenotype}_combined.clumped

  # Create files of SNPs meeting several p-value thresholds. Files numbered 0-4, corresponding values above.
  for threshold in 0 1 2 3 4
  do
    # Further filter clumped SNPs using p-value thresholds (removes multiallelic SNPs)
    python scripts2/06c_filter_snps_for_prs.py \
      data/gwas_results/${phenotype}_combined.clumped \
      --threshold ${thresholds[$threshold]} \
      --output data/prs/${phenotype}_threshold_${threshold}.txt
  done

  # Create combined GWAS result files for each phenotype
  python scripts2/06d_combine_glm_threshold_4.py \
    $pipe/data/gwas_results/${phenotype}.chr*.${phenotype}.glm.linear \
    --keep data/prs/${phenotype}_threshold_4.txt \
    --output data/gwas_results/${phenotype}_combined.glm.linear
done

####################################

# Combine all SNPs meeting the least extreme threshold for any phenotype
cat data/prs/*_threshold_4.txt | uniq > data/prs/all_prs_snps.txt

# Extract any SNPs meeting the least extreme criteria from all chromosomes
rm -f data/prs/prs_merged_list.txt


for chromosome in $(seq 1 22);
do
  $plink2\
    --memory 35000 \
    --bgen /rigel/mfplab/projects/ukb_hakhamanesh/imputed/bgen_files/_001_ukb_imp_chr${chromosome}_v2.bgen ref-first \
    --sample /rigel/mfplab/projects/ukb_hakhamanesh/imputed/bgen_files/ukb_imp.sample \
    --keep data/ukb_population/AFR_all.txt data/ukb_population/AMR_all.txt data/ukb_population/EAS_all.txt data/ukb_population/EUR_test.txt data/ukb_population/SAS_all.txt 
    --make-bed \
    --out data/prs/${chromosome}_temp

  $plink2 \
    --memory 35000 \
    --bfile data/prs/${chromosome}_temp \
    --exclude data/ukb_bbj_pgen/${chromosome}_duplicates.txt \
    --extract data/prs/all_prs_snps.txt \
    --make-bed \
    --out data/prs/${chromosome}_temp1

  rm data/prs/${chromosome}_temp.*

  printf "data/prs/%s_temp1\n" $chromosome >> data/prs/prs_merged_list.txt
done

$plink\
  --merge-list data/prs/prs_merged_list.txt \
  --make-bed \
  --memory 35000 \
  --out data/prs/test_UKBB_imp_BBJ_GWAS_1KG_EUR_LD_merged

rm data/prs/*_temp1.*
