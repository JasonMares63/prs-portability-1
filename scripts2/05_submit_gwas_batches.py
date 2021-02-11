import os
import pathlib
import pandas as pd


def submit_batch_fst_computation(phenotype, temp_dir):
    """
    Mass batch of phenotype
    """

    script_path = temp_dir.joinpath(f'run_{phenotype}.sh')
    script_string = (
        "#!/bin/bash\n"
        "#SBATCH --account=mfplab\n"
        f"#SBATCH --job-name=GWAS_{phenotype}\n"
        "#SBATCH -c 1\n"
        "#SBATCH --time=00:45:00\n"
        "#SBATCH --mem-per-cpu=8gb\n\n"
        "plink2='/rigel/mfplab/users/jm4454/plink/plink2'\n"
        f"phenotype={phenotype}\n"
        f"for chromosome in $(seq 1 22);\n"
        "do\n"
        "  $plink2 "
        "  --bfile data/ukb_filtered/chr${chromosome} "
        "  --keep data/ukb_populations/${phenotype}1.txt "
        "  --pheno data/phenotypes/full_phenotypes.pheno "
        "  --pheno-name $phenotype "
        "  --require-pheno $phenotype "
        "  --covar data/ukb_merged/covar_all_samples.covar "
        "  --covar-name age sex_covar age_sq age_sex age_sq_sex $(printf 'PC%i_AVG ' $(seq 1 20)) "
        "  --require-covar age sex_covar age_sq age_sex age_sq_sex $(printf 'PC%i_AVG ' $(seq 1 20)) "
        "  --vif 100000 "
        "  --memory 35000 "
        "  --glm hide-covar "
        "  --out data/gwas_results/${phenotype}.chr${chromosome}\n "
        "done\n"
        # Delete the script itself
        f"rm {script_path.as_posix()}")
    with open(script_path, 'w') as f:
        f.write(script_string)
    os.system(f'sbatch {script_path.as_posix()}')


def main():
    temp_dir = pathlib.Path('data/gwas_results/')
#    temp_dir.mkdir(exist_ok=True)

    phenotypes = ["BMI","Lymphocyte", "Height", "Eosinophil", "MCH", "MCV", "Monocyte", "Platelet", "RBC", "WBC"]
    for phenotype in phenotypes:
        submit_batch_fst_computation(phenotype, temp_dir)


if __name__ == "__main__":
    main()
