#!/bin/sh
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=classify      # The job name.
#SBATCH -c 5                     # The number of cpu cores to use.
#SBATCH --time=60:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=5gb        # The memory the job will use per cpu core.

set -e

module load anaconda/3-4.2.0
source activate geno

python3 scripts/06b_classify_ukb.py
