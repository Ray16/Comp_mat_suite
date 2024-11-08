import os
from glob import glob
from pathlib import Path
from ase.io import read, write

for traj_path in glob('../MOF_*_CO2/*K/*.traj'):
    print(f'working on {traj_path.split('/')[-1]}')
    parent_folder = Path(traj_path).parent.resolve().absolute()
    sys_name = Path(traj_path).stem
    # Load the trajectory file
    trajectory = read(traj_path, index=':')

    # Write the trajectory to an XYZ file
    write(os.path.join(parent_folder,f'{sys_name}.xyz'), trajectory)