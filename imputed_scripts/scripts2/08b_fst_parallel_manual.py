import os 
import pathlib

temp_dir = pathlib.Path('temp_fst_path/')
for group_label in range(102,121):
    script_path = temp_dir.joinpath(f'run_{group_label}.sh')
    os.system(f'sbatch {script_path.as_posix()}')
