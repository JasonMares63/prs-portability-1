#!/bin/sh
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=merge         # The job name.
#SBATCH -c 6                     # The number of cpu cores to use.
#SBATCH --time=60:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=3gb        # The memory the job will use per cpu core.

/rigel/mfplab/users/mnz2108/plink/plink \
  --merge-list data/merged_list.txt \
  --make-bed \
  --out data/merged/merged

/rigel/mfplab/users/mnz2108/plink/plink2 \
  --bfile data/merged/merged \
  --make-pgen \
  --out data/merged/merged
