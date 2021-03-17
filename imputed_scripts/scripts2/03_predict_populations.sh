#!/bin/sh
#
#SBATCH --account=mfplab
#SBATCH --job-name=classify_populations
#SBATCH -c 3
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=10gb

set -e

module load anaconda/3-4.2.0
source activate prs1

# Train a classifier on 1000 Genomes to predict super population labels
# Save classifier and evaluate it on the UK Biobank PCs
python scripts2/03a_classify_ukb.py

# Create files with sample IDs for each group
# Also does train/test splitting for the "target" individuals.
python scripts2/03b_separate_populations.py

# Rename erroneously named file data/ukb_populations/EUR_all = all training EUR

# Create plots of PCAs of 1000 Genomes population classifier on 1KG and  UK Biobank
Rscript scripts2/03c_plot_pca.R
