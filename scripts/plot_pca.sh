#!/bin/sh
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=pcaPlot       # The job name.
#SBATCH -c 1                     # The number of cpu cores to use.
#SBATCH --time=60:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=3gb        # The memory the job will use per cpu core.

module load anaconda/3-4.2.0
source activate geno

Rscript scripts/_plot_pca.R
