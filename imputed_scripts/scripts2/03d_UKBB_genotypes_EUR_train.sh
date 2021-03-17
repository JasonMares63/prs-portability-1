#!/bin/bash
#
#SBATCH --account=mfplab
#SBATCH --job-name=EUR_train_geno
#SBATCH -c 3
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=8gb


set -e

plink='/rigel/mfplab/users/jm4454/plink/plink'
plink2='/rigel/mfplab/users/jm4454/plink/plink2'


for chromosome in $(seq 1 22);
do
   $plink2 \
      --bfile data/ukb_merged/chr${chromosome} \
      --keep data/ukb_populations/LD_EUR_train.txt \
      --make-bed \
      --out data/LDpred/LD_EUR_train_${chromosome}

   printf "data/LDpred/LD_EUR_train_%s\n" $chromosome >> data/LDpred/LD_EUR_train_merged_list.txt
done

$plink \
  --merge-list data/LDpred/LD_EUR_train_merged_list.txt \
  --make-bed  \
  --out data/LDpred/LD_EUR_merged

#$plink \
#  --bfile data/ukb_merged/merged \
#  --keep data/ukb_populations/EUR_all.txt \
#  --make-bed \
#  --out data/ukb_merged/EUR_all
