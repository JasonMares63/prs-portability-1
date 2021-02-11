#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=merged
#SBATCH -c 10
#SBATCH --time=0-00:30:00
#SBATCH --mem-per-cpu=6gb

set -e

for i in $(seq 1 22);
do
     printf 'data/ukb_bbj_pgen/%s\n'  $i >> data/ukb_bbj_pgen/merged_list.txt
done

plink='/rigel/mfplab/users/jm4454/plink/plink'

$plink \
  --memory 35000 \
  --merge-list data/ukb_bbj_pgen/merged_list.txt \
  --make-bed \
  --out data/ukb_bbj_pgen/merged_fst


