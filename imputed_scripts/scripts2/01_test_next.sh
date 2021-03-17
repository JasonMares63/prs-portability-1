#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=filter_attempt_impute22
#SBATCH -c 1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=10gb

set -e
plink="/rigel/mfplab/users/jm4454/plink/plink"
plink2="/rigel/mfplab/users/jm4454/plink/plink2"

wb='/rigel/mfplab/projects/ukb_hakhamanesh/imputed/plink_neale'
nonwb='/rigel/mfplab/projects/ukb_hakhamanesh/imputed/plink_neale_nonWB'

i=22

$plink2 \
    --pfile data/ukb_filtered/wb_chr$i \
    --extract data/ukb_filtered/wb_chr${i}.prune.in \
    --remove /rigel/mfplab/projects/prs-portability/data/ukb_meta/excluded_samples.sam \
    --make-bed \
    --out data/ukb_filtered/wb_chr$i

rm data/ukb_filtered/wb_chr${i}.pgen
rm data/ukb_filtered/wb_chr${i}.psam
rm data/ukb_filtered/wb_chr${i}.pvar
