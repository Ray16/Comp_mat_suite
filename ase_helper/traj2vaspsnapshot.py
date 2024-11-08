import os
from glob import glob
from pathlib import Path
from ase.io import read, write

traj_path = 'path_to_traj_files'
output_dir = 'path_to_output_files'
# define interval of snapshots
interval = 100

for traj_path in glob(traj_path):
    print(f'working on {traj_path.split('/')[-1]}')
    parent_folder = Path(traj_path).parent.resolve().absolute()
    sys_name = Path(traj_path).stem
    T = sys_name.split('_')[-1]
    sys = traj_path.split('/')[1]
    # Load traj file
    trajectory = read(traj_path, index=f'::{interval}')
    for idx, snapshot in enumerate(trajectory):
        time_ps = 0.01 * 100 * idx
        os.makedirs(output_dir,exist_ok = True)
        # Write the last frame of the trajectory to an poscar file
        write(os.path.join(output_dir,f'{int(time_ps)}ps.vasp'), snapshot, format='vasp')