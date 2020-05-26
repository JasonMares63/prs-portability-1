#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=fst
#SBATCH -c 1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=1gb

set -e

module load anaconda/3-4.2.0
source activate prs

python 08a_make_fst_clusters.py
python scripts/08b_fst_parallel.py

# source deactivate
