import os
import pandas as pd
from glob import glob
from ase import units
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from ase.md import MDLogger
from ase.io import read, write
from ase.optimize import BFGS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution,Stationary
import numpy as np
import time
from time import perf_counter
from mace.calculators import MACECalculator
from ase.constraints import FixAtoms

calc = MACECalculator(model_paths='/project/lgagliardi/ray/Reactive_MLFF/mace_models/2023-12-03-mace-128-L1_epoch-199.model', device='cuda')

structure_dir = '/project/lgagliardi/ray/Raman_spectrum/packmol_sequential_insertion/Zr-MOF-801'
all_structure = []
all_energy = []
for structure in sorted(glob(structure_dir+'/*.cif')):
    all_structure.append('Zr_'+structure.split('/')[-1].split('.cif')[0])
    structure_name = structure.split('/')[-1].split('.cif')[0]
    atoms = read(structure)
    atoms.calc = calc

    slab_e = atoms.get_potential_energy()
    all_energy.append(slab_e)
df = pd.DataFrame({'structure':all_structure,'energy':all_energy})
df.to_csv('structure_energy_Zr.csv',index=False)