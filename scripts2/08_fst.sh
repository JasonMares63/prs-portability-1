#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=fst
#SBATCH -c 1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=1gb

set -e

module load anaconda/3-4.2.0
source activate prs1


#cp data/ukb_merged/merged.fam data/ukb_merged/merged_copy.fam
#python scripts2/08a_make_fst_clusters.py
Rscript scripts2/08b_edit_merged_fam.R
python scripts/08b_fst_parallel.py
#Rscript scripts2/08d_final_fst_formatting.R
# source deactivate
