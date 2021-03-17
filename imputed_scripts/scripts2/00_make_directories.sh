#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=prepare_directories
#SBATCH -c 1
#SBATCH --time=5:00
#SBATCH --mem-per-cpu=1gb

mkdir data/ambiguous_indel_snps \
      data/intersecting_filtered \
      data/kgp_filtered \
      data/kgp_merged \
      data/kgp_meta \
      data/ukb_filtered \
      data/ukb_merged \
      data/ukb_meta \
      data/ukb_populations \
      data/models \
      data/phenotypes \
      data/gwas_results \
      data/prs \
      data/kgp_populations \
      data/fst \ 
      data/LDpred \
      data/LDpred/prs \
      data/prs_comparisons \ 
      img \
      ldpred_scripts
