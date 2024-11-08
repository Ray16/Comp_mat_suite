import os
import pandas as pd
from tqdm import tqdm
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
import matplotlib.pyplot as plt
from ase.io import read, write
import matplotlib.pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import CrystalNN

# The indices are +1 because of indexing from 1 instead of 0

traj_path = '/project/lgagliardi/ray/COF9/C-N-bond-length/traj/liquid_amine_NPT_NVT/traj_npt_298K_1bar_step4.xyz'

sys = traj_path.split('/')[1]

traj_ase = read(traj_path,index=':')
print('Loaded trajectory in ase')

# load traj
adaptor = AseAtomsAdaptor()
traj = []

for atoms in tqdm(traj_ase,total=len(traj_ase)):
    pymatgen_structure = adaptor.get_structure(atoms)
    traj.append(pymatgen_structure)
print('Loaded trajectory for pymatgen')

#### get C_idx in CO2 and N_idx in primary and secondary amines
# read first snapshot
structure = traj[0]
list_symbol = [atom.species_string for atom in structure]

# indices of N
list_N_idx = []
for idx,symbol in enumerate(list_symbol):
    if symbol == 'N':
        list_N_idx.append(idx)

# indices of C
list_C_idx = []
for idx,symbol in enumerate(list_symbol):
    if symbol == 'C':
        list_C_idx.append(idx)

cutoff=2.8
nn = CrystalNN(search_cutoff=cutoff)

# detect C atoms in CO2
CO2_C_idx = []
CO2_C_label = []
for C_idx in list_C_idx:
    nn_info = nn.get_nn_info(structure,C_idx)
    if len(nn_info) == 2:
        list_site_symbol = [site_dict['site'].species_string for site_dict in nn_info]
        # CO2
        if list_site_symbol.count('O') == 2:
            CO2_C_idx.append(C_idx)
            CO2_C_label.append(structure[C_idx].label)

print(f'Number of CO2 molecules: {len(CO2_C_idx)}')

# detect N atoms in primiary amines and secondary amines
primary_amine_N_idx = []
primary_amine_N_label = []
secondary_amine_N_idx = []
secondary_amine_N_label = []
for N_idx in list_N_idx:
    nn_info = nn.get_nn_info(structure,N_idx)
    if len(nn_info) == 3:
        list_site_symbol = [site_dict['site'].species_string for site_dict in nn_info]
        # primary amine
        if list_site_symbol.count('C') == 1 and list_site_symbol.count('H') == 2:
            primary_amine_N_idx.append(N_idx)
            primary_amine_N_label.append(structure[N_idx].label)
        # secondary amine
        if list_site_symbol.count('C') == 2 and list_site_symbol.count('H') == 1:
            secondary_amine_N_idx.append(N_idx)
            secondary_amine_N_label.append(structure[N_idx].label)

print(f'Number of primary amines: {len(primary_amine_N_idx)}')
print(f'Number of secondary amines: {len(secondary_amine_N_idx)}')

os.makedirs(f'atom_idx/{sys}',exist_ok=True)

# output C_idx for C in CO2 and N_idx for N in amines
dict_idx = {'C_CO2_idx':CO2_C_idx,'C_CO2_label':CO2_C_label,'N_primary_amine_idx':primary_amine_N_idx,'N_primary_amine_label':primary_amine_N_label,'N_secondary_amine_idx':secondary_amine_N_idx,'N_secondary_amine_label':secondary_amine_N_label}
with open(f'atom_idx/{sys}/C_N_idx_label.pkl','wb') as f:
    pickle.dump(dict_idx,f)
print('Saved C and N indices and labels to pickle file')

#### output and plot C-N distance ####
os.makedirs(f'output/{sys}',exist_ok=True)
os.makedirs(f'plot/{sys}',exist_ok=True)

# calculate C-N distance
for dict_C_idx, C_idx in tqdm(enumerate(dict_idx['C_CO2_idx'])):
    plt.figure()
    min_list_C_N_distance = []
    min_list_C_N_distance_N_idx = []
    # iterature through the snapshots
    for atoms in traj:
        # iterature through the primary amine and secondary amine N indices
        N_indices = dict_idx['N_primary_amine_idx']+dict_idx['N_secondary_amine_idx']
        list_C_N_distance_snapshot = []
        for dict_N_idx,N_idx in enumerate(N_indices):
            # calculate C-N distance
            C_N_distance = atoms.lattice.get_distance_and_image(atoms[N_idx].frac_coords, atoms[C_idx].frac_coords)[0]
            list_C_N_distance_snapshot.append(C_N_distance)
        min_list_C_N_distance.append(min(list_C_N_distance_snapshot))
        min_list_C_N_distance_N_idx.append(N_indices[list_C_N_distance_snapshot.index(min(list_C_N_distance_snapshot))])
    # save C-N distance to file
    df = pd.DataFrame({'Timestep':range(len(traj)),'min_C_N_distance':min_list_C_N_distance,'min_C_N_distance_N_idx':min_list_C_N_distance_N_idx})
    df.to_csv(f'output/{sys}/C{C_idx}.csv',index=False)

    plt.plot(range(len(traj)),min_list_C_N_distance)
    plt.axhline(1.4,linestyle='--',linewidth=2,c='gray')
    plt.tight_layout()
    plt.xlabel('Timestep')
    plt.ylabel('Minimum C-N distance (Å)')
    plt.savefig(f"plot/{sys}/C{C_idx}.jpg",dpi=300,bbox_inches='tight')
    plt.close()
