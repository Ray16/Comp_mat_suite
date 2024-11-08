import argparse
from glob import glob
from ase.io import read, write

parser = argparse.ArgumentParser()
parser.add_argument('--idx',type=str)
args = parser.parse_args()

cif_file = f'h2o_{args.idx}.cif'
cif_name = cif_file.split('.cif')[0]
atoms = read(cif_file)
write(cif_name+'.xyz', atoms)