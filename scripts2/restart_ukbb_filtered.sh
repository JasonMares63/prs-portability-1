#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=make_filtered_genotypes_files
#SBATCH -c 5
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=10gb

set -e

plink="/rigel/mfplab/users/jm4454/plink/plink"
plink2="/rigel/mfplab/users/jm4454/plink/plink2"
module load anaconda/3-4.2.0
source activate prs1

for i in $(seq 1 22);
do
  # Identify indels and ambiguous variants and write them to a file
  python scripts2/01a_get_ambiguous_indel_snps.py \
    /rigel/mfplab/projects/ukb_hakhamanesh/genotyped/chr${i}.bim \
    -o data/ambiguous_indel_snps/chr${i}.snps

  # Remove the ambiguous or indel variants and write the resulting SNPs to a file
  python scripts2/01b_remove_ambiguous_indel_snps.py \
    /rigel/mfplab/projects/ukb_hakhamanesh/genotyped/chr${i}.bim \
    -r data/ambiguous_indel_snps/chr${i}.snps \
    -o data/intersecting_filtered/chr${i}.snps

  # Convert from Plink 1 to Plink 2 and compute the SNPs to be dropped
  $plink2 \
    --bfile  /rigel/mfplab/projects/ukb_hakhamanesh/genotyped/chr${i}  \
    --maf 0.05 \
    --geno 0.01 \
    --indep-pairwise 1000 kb 1 0.2 \
    --extract data/intersecting_filtered/chr${i}.snps \
    --make-pgen \
    --out data/ukb_filtered/chr$i

  # Write a new file using only filtered SNPs for 1000 Genomes (Plink 1 format)
  $plink2 \
    --pfile data/ukb_filtered/chr$i \
    --extract data/ukb_filtered/chr${i}.prune.in \
    --make-bed \
    --out data/ukb_filtered/chr$i

  # Append the output file path to a new file (for merging them all below)
  printf "data/ukb_filtered/chr%s\n" $i >> data/ukb_merged/ukb_merged_list.txt
done


$plink \
  --merge-list data/ukb_merged/ukb_merged_list.txt \
  --make-bed \
  --out data/ukb_merged/merged

$plink2 \
  --bfile data/ukb_merged/merged \
  --make-pgen \
  --out data/ukb_merged/merged
