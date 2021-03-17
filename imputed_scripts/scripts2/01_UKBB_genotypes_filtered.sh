#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=make_filtered_genotypes_files
#SBATCH -c 5
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=16gb

set -e

plink="/rigel/mfplab/users/jm4454/plink/plink"
plink2="/rigel/mfplab/users/jm4454/plink/plink2"

wb='/rigel/mfplab/projects/ukb_hakhamanesh/imputed/plink_neale'
nonwb='/rigel/mfplab/projects/ukb_hakhamanesh/imputed/plink_neale_nonWB'

module load anaconda
source activate prs1

for i in $(seq 1 23);
do
  # Identify indels and ambiguous variants and write them to a file
  python scripts2/01a_get_ambiguous_indel_snps.py \
    ${wb}/chr${i}.bim \
    -o data/ambiguous_indel_snps/chr${i}.snps

  # Convert from Plink 1 to Plink 2 and compute the SNPs to be dropped
  $plink2 \
    --bfile ${wb}/chr${i}  \
    --keep ${wb}/ukb_white_indivs_neal.fam \
    --remove /rigel/mfplab/projects/prs-portability/data/ukb_meta/excluded_samples.sam \
    --exclude data/ambiguous_indel_snps/chr${i}.snps \
    --maf 0.01 \
    --geno 0.01 \
    --indep-pairwise 1000 kb 1 0.2 \
    --make-pgen \
    --out data/ukb_filtered/wb_chr$i

  # Write a new file using only filtered SNPs for 1000 Genomes (Plink 1 format)
  $plink2 \
    --pfile data/ukb_filtered/wb_chr${i} \
    --extract data/ukb_filtered/wb_chr${i}.prune.in \
    --make-bed \
    --out data/ukb_filtered/wb_chr${i}

  rm data/ukb_filtered/wb_chr${i}.pgen
  rm data/ukb_filtered/wb_chr${i}.pvar
  rm data/ukb_filtered/wb_chr${i}.psam 

  $plink2 \
    --bfile ${nonwb}/chr${i} \
    --keep ${nonwb}/ukb_filter_indivs_nowhiteB_unrelated.txt \
    --remove /rigel/mfplab/projects/prs-portability/data/ukb_meta/excluded_samples.sam \
    --extract data/ukb_filtered/wb_chr${i}.prune.in \
    --make-bed \
    --out data/ukb_filtered/nonwb_chr${i}


  # Append the output file path to a new file (for merging them all below)
  printf "data/ukb_filtered/wb_chr%s\n" $i >> data/ukb_merged/ukb_merged_list${i}.txt
  printf "data/ukb_filtered/nonwb_chr%s\n" $i >> data/ukb_merged/ukb_merged_list${i}.txt

  $plink \
    --merge-list data/ukb_merged/ukb_merged_list${i}.txt \
    --make-bed \
    --out data/ukb_merged/chr${i}

  rm data/ukb_filtered/wb_chr${i}.b*
  rm data/ukb_filtered/nonwb_chr${i}.b*
  rm data/ukb_filtered/*.fam
  rm data/ukb_filtered/*.log

  printf "data/ukb_merged/chr$%s\n" $i >> data/ukb_merged/ukb_merged_list.txt
done


#$plink \
#  --merge-list data/ukb_merged/ukb_merged_list.txt \
#  --make-bed \
#  --out data/ukb_merged/merged

#$plink2 \
#  --bfile data/ukb_merged/merged \
#  --make-pgen \
#  --out data/ukb_merged/merged


