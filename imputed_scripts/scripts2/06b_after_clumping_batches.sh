#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=post-clump
#SBATCH -c 1
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=8gb

set -e

module load anaconda/3-4.2.0
source activate prs1

plink="/rigel/mfplab/users/jm4454/plink/plink"

cat data/prs/*_threshold_4.txt | uniq > data/prs/all_prs_snps.txt

for chromosome in $(seq 1 22);
do
  $plink \
    --bfile data/ukb_merged/chr${chromosome} \
    --extract data/prs/all_prs_snps.txt \
    --memory 35000 \
    --make-bed \
    --out data/prs/chr${chromosome}_temp

  printf "data/prs/chr%s_temp\n" $chromosome >> data/prs/prs_merged_list.txt
done

# Combine extracted SNPs across chromosomes into a single Plink 1 file of all
#  SNPs that meet the least extreme level of significance for any of the traits.
$plink \
  --merge-list data/prs/prs_merged_list.txt \
  --make-bed \
  --memory 35000 \
  --out data/prs/CT_merged

rm data/prs/chr*_temp.*
