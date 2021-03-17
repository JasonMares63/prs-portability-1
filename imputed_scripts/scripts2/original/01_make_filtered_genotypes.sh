#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=make_filtered_genotypes_files
#SBATCH -c 5
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=10gb

# Fail if any command fails
set -e

module load anaconda/3-4.2.0
source activate prs

# 1000 Genomes Project
###############################################################################
for i in $(seq 1 22);
do
  # Identify indels and ambiguous variants and write them to a file
  python scripts/01a_get_ambiguous_indel_snps.py \
    /rigel/mfplab/users/hsm2137/1000_genome/data/plink_files/chr${i}.bim \
    -o data/ambiguous_indel_snps/chr${i}.snps

  # Remove the ambiguous or indel variants and write the resulting SNPs to a file
  python scripts/01b_remove_ambiguous_indel_snps.py \
    /rigel/mfplab/users/hsm2137/1000_genome/data/plink_files/ukb/chr${i}_shared.snps \
    -r data/ambiguous_indel_snps/chr${i}.snps \
    -o data/intersecting_filtered/chr${i}.snps

  # Convert from Plink 1 to Plink 2 and compute the SNPs to be dropped
  /rigel/mfplab/users/mnz2108/plink/plink2 \
    --bfile /rigel/mfplab/users/hsm2137/1000_genome/data/plink_files/chr$i \
    --maf 0.05 \
    --geno 0.01 \
    --indep-pairwise 1000 kb 1 0.2 \
    --extract data/intersecting_filtered/chr${i}.snps \
    --make-pgen \
    --out data/kgp_filtered/chr$i

  # Write a new file using only filtered SNPs for 1000 Genomes (Plink 1 format)
  /rigel/mfplab/users/mnz2108/plink/plink2 \
    --pfile data/kgp_filtered/chr$i \
    --extract data/kgp_filtered/chr${i}.prune.in \
    --make-bed \
    --out data/kgp_filtered/chr$i

  # Append the output file path to a new file (for merging them all below)
  printf "data/kgp_filtered/chr%s\n" $i >> data/kgp_merged/kgp_merged_list.txt
done

# Merge all the 1000 Genomes files into a single Plink 1 file
# Currently, Plink 2 does not appear to support this merge operation.
/rigel/mfplab/users/mnz2108/plink/plink \
  --merge-list data/kgp_merged/kgp_merged_list.txt \
  --make-bed \
  --out data/kgp_merged/merged

# Convert the merged Plink 1 file to the Plink 2 format (pgen/pvar/psam)
/rigel/mfplab/users/mnz2108/plink/plink2 \
  --bfile data/kgp_merged/merged \
  --make-pgen \
  --out data/kgp_merged/merged
###############################################################################


# UK Biobank
###############################################################################
for i in $(seq 1 22);
do
  # Make a Plink 1 representation of the chromosome using only the filtered intersecting SNPs
  /rigel/mfplab/users/mnz2108/plink/plink2 \
    --memory 40000 \
    --bgen /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/_001_ukb_imp_chr${i}_v2.bgen ref-first \
    --sample /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/ukb_imp.sample \
    --extract data/kgp_filtered/chr${i}.prune.in \
    --make-bed \
    --out data/ukb_filtered/chr${i}

  # Add the file name to a running list (for easier merging below)
  printf "data/ukb_filtered/chr%s\n" $i >> data/ukb_merged/ukb_merged_list.txt
done

# Merge the Plink 1 format files into one combined Plink 1 file
/rigel/mfplab/users/mnz2108/plink/plink \
  --memory 40000 \
  --merge-list data/ukb_merged/ukb_merged_list.txt \
  --remove data/ukb_meta/excluded_samples.sam \
  --make-bed \
  --out data/ukb_merged/merged

# Convert the merged Plink 1 file to Plink 2 format (.pgen/.pvar, etc.)
/rigel/mfplab/users/mnz2108/plink/plink2 \
  --bfile data/ukb_merged/merged \
  --make-pgen \
  --out data/ukb_merged/merged
###############################################################################
