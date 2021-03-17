#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=make_filtered_1kg
#SBATCH -c 5
#SBATCH --time=01-00:00:00
#SBATCH --mem-per-cpu=10gb

set -e
#wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz 
#wget https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples.20130502.ALL.ped 
#mv integrated_call_samples.20130502.ALL.ped data/kgp_meta
#mv ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz data/kgp_meta
#gzip -d data/kgp_meta/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'


for i in $(seq 1 22);
do
  # Using existing 1KG file which includes all MAF 1%, genotype 1%, and pruned intersecting SNPs with UKBB.
  # From existing variants, extract those that were filtered from UKBB genotype data.
  $plink2 \
  --pfile /rigel/mfplab/projects/prs-portability/data/kgp_filtered/chr$i \
  --extract data/ukb_filtered/wb_chr${i}.prune.in \
  --make-bed \
  --out data/kgp_filtered/chr${i}

  # Add the file name to a running list (for easier merging below)
  printf "data/kgp_filtered/chr%s\n" $i >> data/kgp_merged/kgp_merged_list.txt

done

$plink \
  --merge-list data/kgp_merged/kgp_merged_list.txt \
  --make-bed \
  --out data/kgp_merged/merged

rm data/kgp_filtered/chr*

#$plink2 \
#  --bfile data/kgp_merged/merged \
#  --make-pgen \
#  --out data/kgp_merged/merged
