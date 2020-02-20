#!/bin/bash
#
#SBATCH --account=mfplab         # The account name for the job.
#SBATCH --job-name=separate      # The job name.
#SBATCH -c 8                     # The number of cpu cores to use.
#SBATCH --time=60:00             # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=5gb        # The memory the job will use per cpu core.

set -e

module load anaconda/3-4.2.0
source activate ../envs

# Create files with sample IDs for each group
# Also does train/test splitting for the "target" individuals.
python scripts/07b_separate_populations.py
source deactivate

# Create new pgen files for each population
for population in AFR AMR EAS EUR SAS
do
echo $population
/rigel/mfplab/users/mnz2108/plink/plink2 --pfile data/ukb_filtered/merged --keep data/ukb_populations/${population}_all.txt  --make-pgen --out data/ukb_populations/${population}_all
done
