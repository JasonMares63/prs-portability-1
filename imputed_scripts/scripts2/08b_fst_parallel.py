import itertools
import os
import pathlib

import pandas as pd


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def submit_batch_fst_computation(group_label, test_fids, temp_dir):
    """
    Batch Fst computations. We have ~32,000 test individuals and would rather
    submit O(100) jobs than 32,000 jobs.
    """
    formatted_fids = ' '.join([str(i) for i in test_fids])
    script_path = temp_dir.joinpath(f'run_{group_label}.sh')
    script_string = (
        "#!/bin/bash\n"
        "#SBATCH --account=mfplab\n"
        f"#SBATCH --job-name=fst{group_label}\n"
        "#SBATCH -c 1\n"
        "#SBATCH --time=5-00:00:00\n"
        "#SBATCH --mem-per-cpu=8gb\n\n"
        f"group={group_label}\n"
        f"for test_fid in {formatted_fids}\n"
        "do\n"
        "  /rigel/mfplab/users/jm4454/plink/plink "
        "  --bim data/ukb_merged/merged.bim "
        "  --bed data/ukb_merged/merged.bed "
        "  --fam data/ukb_merged/merged_fst.fam "
        "  --family "
        "  --fst "
        "  --keep-cluster-names train ${test_fid} "
        "  --out data/fst/fst${test_fid}\n\n"
        "  rm data/fst/fst${test_fid}.fst data/fst/fst${test_fid}.nosex\n"
        "  echo ${test_fid} >> data/fst/fst${group}.est\n"
        "  grep 'Fst estimate:' data/fst/fst${test_fid}.log >> data/fst/fst${group}.est\n"
        "  rm data/fst/fst${test_fid}.log\n"
        "done\n"
        # Delete the script itself
        f"rm {script_path.as_posix()}")
    with open(script_path, 'w') as f:
        f.write(script_string)
#    os.system(f'sbatch {script_path.as_posix()}')


def main():
    temp_dir = pathlib.Path('temp_fst_path/')
    temp_dir.mkdir(exist_ok=True)

    all_test = (
        pd.read_csv('data/ukb_merged/merged_fst.fam', sep='\t')
        .loc[lambda df: df['#FID'] != 'train', '#FID']
        .values
        .tolist()
    )

    # all_test = all_test[:5]

    group_size = 40
    for i, test_fid_group in enumerate(grouper(all_test, group_size, '')):
        submit_batch_fst_computation(i, test_fid_group, temp_dir)


if __name__ == "__main__":
    main()
