import os
import pathlib
import pandas as pd


def submit_batch_fst_computation(phenotype, temp_dir):
    """
    Mass batch of clumping across phenotype
    """

    script_path = temp_dir.joinpath(f'run_{phenotype}_clump.sh')
    script_string = (
        "#!/bin/bash\n"
        "#SBATCH --account=mfplab\n"
        f"#SBATCH --job-name=CT_{phenotype}\n"
        "#SBATCH -c 1\n"
        "#SBATCH --time=01:00:00\n"
        "#SBATCH --mem-per-cpu=8gb\n\n"
	"set -e \n"
        "plink='/rigel/mfplab/users/jm4454/plink/plink' \n"
	"pipe='/rigel/mfplab/projects/prs-portability' \n"
        f"phenotype={phenotype}\n"
	"thresholds=(5e-8 1e-6 1e-4 1e-3 1e-2) \n"
  	"# Clump GWAS results \n"
  	"for chromosome in $(seq 1 22); \n"
  	"do \n"
   	"# Convert the GWAS output file from the Plink 2 to Plink 1 \n"
        "#python scripts2/06b_convert_plink2_glm_to_plink1.py "
      	"#data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.linear "
      	"#--output data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc \n"

    	"# Clump GWAS results using LD from 1000 Genomes EUR super-population \n"
    	"$plink "
      	"--memory 35000 "
      	"--bfile data/kgp_populations/EUR_${chromosome} "
      	"--clump $pipe/data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.glm.assoc "
      	"--clump-p1 0.01 "
      	"--clump-p2 1 "
      	"--clump-r2 0.5 "
      	"--clump-kb 250 "
      	"--out data/gwas_results/${phenotype}.chr${chromosome}.${phenotype} \n"

    	"# `--clump` automatically outputs plink.* files. Move and rename these \n"
	"#    mv plink.clumped data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.clumped \n"
	"#    mv plink.log data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.clumped.log \n"
	"#    mv plink.nosex data/gwas_results/${phenotype}.chr${chromosome}.${phenotype}.clumped.nosex \n"
  	"done \n"

  	"# Combine clumped SNPs across chromosomes \n"
  	"head -n 1 data/gwas_results/${phenotype}.chr1.${phenotype}.clumped > "
    	"data/gwas_results/${phenotype}_combined.clumped \n"
  	"tail -n +2 -q data/gwas_results/${phenotype}.chr*.${phenotype}.clumped >> " 
    	"data/gwas_results/${phenotype}_combined.clumped \n"

  	"# Create files of SNPs meeting several p-value thresholds. Files numbered 0-4 \n"
  	"for threshold in 0 1 2 3 4 \n"
  	"do \n"
    	"# Further filter clumped SNPs using p-value thresholds (removes multiallelic SNPs) \n"
    	"python scripts2/06c_filter_snps_for_prs.py "
      	"data/gwas_results/${phenotype}_combined.clumped "
      	"--threshold ${thresholds[$threshold]} "
      	"--output data/prs/${phenotype}_threshold_${threshold}.txt \n"
  	"done \n"

  	"# Create combined GWAS result files for each phenotype \n"
  	"python scripts2/06d_combine_glm_threshold_4.py "
    	"$pipe/data/gwas_results/${phenotype}.chr*.${phenotype}.glm.linear "
    	"--keep data/prs/${phenotype}_threshold_4.txt "
    	"--output data/gwas_results/${phenotype}_combined.glm.linear \n "
        
	"# Delete the script itself \n"
        f"rm {script_path.as_posix()} \n"
	)
    with open(script_path, 'w') as f:
        f.write(script_string)
    os.system(f'sbatch {script_path.as_posix()}')


def main():
    temp_dir = pathlib.Path('')
#    temp_dir.mkdir(exist_ok=True)

    phenotypes = ["BMI","Lymphocyte", "Height", "Eosinophil", "MCH", "MCV", "Monocyte", "Platelet", "RBC","WBC"]
    for phenotype in phenotypes:
        submit_batch_fst_computation(phenotype,temp_dir)


if __name__ == "__main__":
    main()


