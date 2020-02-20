#!/bin/bash
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=FilterSNPs    # The job name.
#SBATCH -c 8                    # The number of cpu cores to use.
#SBATCH --time=60:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=3gb        # The memory the job will use per cpu core.

# Fail if any command fails
set -e

module load anaconda/3-4.2.0
source activate ../envs

# cp /rigel/mfplab/users/hsm2137/1000_genome/data/plink_files/chr*.* data/kgp_raw/
# cp /rigel/mfplab/users/hsm2137/1000_genome/data/plink_files/ukb/*_shared.snps data/intersecting_raw/

for i in $(seq 1 22);
do
  python scripts/get_ambiguous_indel_snps.py data/kgp_raw/chr${i}.bim -o data/ambiguous_indel_snps/chr${i}.snps
  python scripts/remove_ambiguous_indel_snps.py data/intersecting_raw/chr${i}_shared.snps -r data/ambiguous_indel_snps/chr${i}.snps -o data/intersecting_filtered/chr${i}.snps

  /rigel/mfplab/users/mnz2108/plink/plink2 \
    --bfile data/kgp_raw/chr$i \
    --maf 0.05 \
    --geno 0.01 \
    --indep-pairwise 1000 kb 1 0.2 \
    --extract data/intersecting_filtered/chr${i}.snps \
    --make-pgen \
    --out data/kgp_filtered/chr$i
done
