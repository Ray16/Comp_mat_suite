import os
import shutil
from glob import glob

charge_mol_calc_dir = 'chargemol_calc'

files_to_copy = ['AECCAR0', 'AECCAR2', 'CHGCAR', 'POTCAR']

for sys in glob('*-MOF-*'):
    target_dir = charge_mol_calc_dir+'/'+sys
    os.makedirs(target_dir,exist_ok=True)
    for file in files_to_copy:
        shutil.copy(os.path.join(sys,file),target_dir)
    # copy run_chargemol_parallel.job and job_control.txt to the corresponding folders
    shutil.copy(os.path.join(charge_mol_calc_dir,'run_chargemol_parallel.job'),target_dir)
    shutil.copy(os.path.join(charge_mol_calc_dir,'job_control.txt'),target_dir)
