#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=prep-PRS
#SBATCH -c 10
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=6gb

set -e
plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'

# Used filtered UK-BB genotypes +  UKBB EUR-clumped SNPs + BBJ-centered GWAS 
for chromosome in $(seq 1 22);
do
  $plink2\
    --pfile /rigel/mfplab/projects/prs-portability/data/ukb_merged/merged \
    --chr ${chromosome} \
    --extract /rigel/mfplab/projects/prs-portability/data/prs/all_prs_snps.txt \
    --exclude /rigel/mfplab/projects/prs-portability/data/ukb_bbj_pgen/${chromosome}_duplicates.txt \
    --memory 35000 \
    --make-bed \
    --out data/prs/chr${chromosome}_temp

#  printf "data/prs/chr%s_temp\n" $chromosome >> data/prs/prs_merged_list.txt
done

# Combine extracted SNPs across chromosomes into a single Plink 1 file of all
#  SNPs that meet the least extreme level of significance for any of the traits.
$plink\
  --merge-list data/prs/prs_merged_list.txt \
  --make-bed \
  --memory 35000 \
  --out data/prs/merged

rm data/prs/chr*_temp.*

$plink2 \
  --bfile data/prs/merged \
  --make-pgen \
  --memory 35000 \
  --out data/prs/merged
