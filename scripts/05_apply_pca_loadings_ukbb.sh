#!/bin/sh
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=pcaApply      # The job name.
#SBATCH -c 10                    # The number of cpu cores to use.
#SBATCH --time=24:00:00          # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=8gb        # The memory the job will use per cpu core.


# Fail if any command fails
set -e

# Copy UK Biobank genotypes to the local directory (to avoid any chance of accidentally messing something up)
cp /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/_001_ukb_imp_chr*_v2.bgen data/ukb_raw/

for i in $(seq 13 22);
do
  # Make a plink 1 representation of the chromosome using only the filtered intersecting SNPs
  /rigel/mfplab/users/mnz2108/plink/plink2 \
    --bgen /rigel/mfplab/users/hsm2137/ukbiobank/data/imputed/bgen_files/_001_ukb_imp_chr${i}_v2.bgen ref-first \
    --sample data/ukb_raw/ukb_imp.sample \
    --extract data/kgp_filtered/chr${i}.prune.in \
    --make-bed \
    --out data/ukb_filtered/chr${i}

  # Add the file name to a running list (for easier merging below)
  printf "data/ukb_filtered/chr%s\n" $i >> data/ukb_merged_list.txt
done
 
# Merge the plink 1 format files into one combined plink 1 file
/rigel/mfplab/users/mnz2108/plink/plink \
  --merge-list data/ukb_merged_list.txt \
  --make-bed \
  --out data/ukb_filtered/merged
 
# Convert the merged plink 1 file to plink 2 format (.pgen/.pvar, etc.)
/rigel/mfplab/users/mnz2108/plink/plink2 \
  --bfile data/ukb_filtered/merged \
  --make-pgen \
  --out data/ukb_filtered/merged
 
# Apply loadings from PCA computed on 1000 Genomes to the merged plink 2 file
/rigel/mfplab/users/mnz2108/plink/plink2 \
  --pfile data/ukb_filtered/merged \
  --read-freq data/merged/merged.acount \
  --score data/merged/merged.eigenvec.allele 2 5 header-read no-mean-imputation \
  --score-col-nums 6-25 \
  --out data/ukb_filtered/projection

#  Apply loadings from PCA on 1000 Genomes to 1000 Genomes (for consistent method)
/rigel/mfplab/users/mnz2108/plink/plink2 \
  --pfile data/merged/merged \
  --read-freq data/merged/merged.acount \
  --score data/merged/merged.eigenvec.allele 2 5 header-read no-mean-imputation \
  --score-col-nums 6-25 \
  --out data/merged/projection

