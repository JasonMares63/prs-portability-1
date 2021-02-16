#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=prep-PRS-ukb-bbj-gwas-ukb-ld
#SBATCH -c 10
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=4gb

set -e

#Central directory of plink not up to-date
plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'
pipe='/rigel/mfplab/projects/prs-portability'
module load anaconda/3-4.2.0
source activate prs1

cat data/prs/*_threshold_4.txt | uniq > data/prs/all_prs_snps.txt

# Extract any SNPs meeting the least extreme criteria from all chromosomes
rm -f data/prs/prs_merged_list.txt


for chromosome in $(seq 1 22);
do
  $plink2 \
    --memory 35000 \
    --bgen /rigel/mfplab/projects/ukb_hakhamanesh/imputed/bgen_files/_001_ukb_imp_chr${chromosome}_v2.bgen ref-first \
    --sample /rigel/mfplab/projects/ukb_hakhamanesh/imputed/bgen_files/ukb_imp.sample \
    --keep data/ukb_populations/AFR_all.txt data/ukb_populations/AMR_all.txt data/ukb_populations/EAS_all.txt data/ukb_populations/EUR_test.txt data/ukb_populations/SAS_all.txt \
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

# Combine extracted SNPs across chromosomes into a single Plink 1 file of all
#  SNPs that meet the least extreme level of significance for any of the traits.
$plink\
  --merge-list data/prs/prs_merged_list.txt \
  --make-bed \
  --memory 35000 \
  --out data/prs/test_UKBB_imp_BBJ_GWAS_1KG_EUR_LD_merged

rm data/prs/*_temp1.*
