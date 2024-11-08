import os
import shutil
import argparse
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

parser = argparse.ArgumentParser()
parser.add_argument('--idx',type=str)
args = parser.parse_args()

structure_dir = '/project/lgagliardi/ray/Raman_spectrum/packmol_sequential_insertion/Zr-MOF-801/'

# backup the original cif file
shutil.copy(f'h2o_{args.idx}.cif', f'h2o_{args.idx}_unoptimized.cif')

structure_cif_name = os.path.join(structure_dir, f'h2o_{args.idx}.cif')
if structure_cif_name not in os.listdir('.'):
    print(f'Optimizing {structure_cif_name}')
    atoms = read(structure_cif_name)
    atoms.calc = calc
    atoms.set_pbc([True, True, True])

    opt = BFGS(atoms)
    opt.run(fmax=0.03, steps=1000)
    if opt.nsteps <= 1000:
        write(structure_cif_name,atoms)