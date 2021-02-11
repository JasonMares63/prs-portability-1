#!/bin/sh
#
#SBATCH --account=mfplab
#SBATCH --job-name=classify_populations
#SBATCH -c 1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=10gb

set -e

module load anaconda/3-4.2.0
source activate prs

# Train a classifier on 1000 Genomes to predict super population labels
# Save classifier and evaluate it on the UK Biobank PCs
python scripts/03a_classify_ukb.py

# Create files with sample IDs for each group
# Also does train/test splitting for the "target" individuals.
python scripts/03b_separate_populations.py

# Create plots of PCAs of 1000 Genomes population classifier on 1kg and  UK Biobank
Rscript scripts/03c_plot_pca.R
