#!/bin/bash
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=prepareMerge  # The job name.
#SBATCH -c 8                     # The number of cpu cores to use.
#SBATCH --time=60:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=3gb        # The memory the job will use per cpu core.

# Fail if any command fails
set -e

for i in $(seq 1 22); 
do
  /rigel/mfplab/users/mnz2108/plink/plink2 \
    --pfile data/kgp_filtered/chr$i \
    --extract data/kgp_filtered/chr${i}.prune.in \
    --make-bed \
    --out data/merged/chr$i
done
