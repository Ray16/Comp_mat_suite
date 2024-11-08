import argparse
from glob import glob
from ase import io
from pymatgen.core import Structure, Lattice
from pymatgen.io.cif import CifWriter

parser = argparse.ArgumentParser()
parser.add_argument('--idx',type=str)
args = parser.parse_args()

#all_xyz_files = glob('*.xyz')

#for xyz_file in all_xyz_files:
xyz_file = f'h2o_{args.idx}.xyz'
MOF_name = xyz_file.split('.xyz')[0]
lattice = Lattice.from_parameters(a=17.95246895, b=17.95246895, c=17.95246895, alpha=90, beta=90, gamma=90)
atoms = io.read(xyz_file)
species = [atom.symbol for atom in atoms]
cartesian_coords = atoms.positions
fractional_coords = lattice.get_fractional_coords(cartesian_coords)
structure = Structure(lattice, species, fractional_coords)

# Write to CIF file
cif_file = f'{MOF_name}.cif'
cif_writer = CifWriter(structure)
cif_writer.write_file(cif_file)