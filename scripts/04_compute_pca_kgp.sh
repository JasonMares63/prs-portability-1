#!/bin/sh
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=pcaCompute    # The job name.
#SBATCH -c 8                     # The number of cpu cores to use.
#SBATCH --time=60:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=5gb        # The memory the job will use per cpu core.

/rigel/mfplab/users/mnz2108/plink/plink2 \
  --pfile data/merged/merged \
  --freq counts \
  --pca allele-wts 20 \
  --out data/merged/merged
