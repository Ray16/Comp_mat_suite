import os
import numpy as np
from glob import glob

def xyz2cif(xyz_file,sys):
    # Extract unit cell information from the XYZ file metadata
    metadata_line = xyz_file[1]
    unit_cell_info = metadata_line.split("unitcell [")[1].split("]")[0]

    # Parse unit cell parameters correctly by removing any braces and splitting by spaces and commas
    unit_cell_params = []

    for part in unit_cell_info.split('},'):
        part = part.replace('{', '').replace('}', '').strip()
        unit_cell_params.append([float(x) for x in part.split()])

    unit_cell_params = np.array(unit_cell_params)

    # Calculate lattice constants
    a = np.linalg.norm(unit_cell_params[0])
    b = np.linalg.norm(unit_cell_params[1])
    c = np.linalg.norm(unit_cell_params[2])
    volume = np.dot(unit_cell_params[0], np.cross(unit_cell_params[1], unit_cell_params[2]))

    # Define the CIF file header and format with an additional charges column
    cif_header_with_charges = """data_generated
    _symmetry_space_group_name_H-M    'P 1'
    _cell_length_a    {:.6f}
    _cell_length_b    {:.6f}
    _cell_length_c    {:.6f}
    _cell_angle_alpha    90.000
    _cell_angle_beta     90.000
    _cell_angle_gamma    90.000
    _symmetry_Int_Tables_number    1
    _chemical_formula_sum  '{}'
    _cell_volume   {:.6f}
    _cell_formula_units_Z    1
    loop_
    _symmetry_equiv_pos_as_xyz
      'x, y, z'
    loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_U_iso_or_equiv
    _atom_site_charge
    """

    # Extract atomic positions, species, and charges, skipping invalid lines
    atom_lines = xyz_file[2:]
    atomic_positions = []
    atomic_species = []
    atomic_charges = []

    for line in atom_lines:
        parts = line.split()
        if len(parts) < 5:  # Ensure that the line has at least 5 elements (element symbol, x, y, z, charge)
            continue
        element = parts[0]
        try:
            x, y, z, charge = map(float, parts[1:])
            atomic_positions.append([x, y, z])
            atomic_species.append(element)
            atomic_charges.append(charge)
        except ValueError:
            # Skip the line if it contains non-numeric data in the coordinates or charge fields
            continue

    # Create CIF file content with charges
    cif_content_with_charges = cif_header_with_charges.format(a, b, c, " ".join(set(atomic_species)), volume)

    # Add atomic positions and charges to CIF
    for i, (element, (x, y, z), charge) in enumerate(zip(atomic_species, atomic_positions, atomic_charges)):
        cif_content_with_charges += f"{element}{i+1} {element} {x/a:.6f} {y/b:.6f} {z/c:.6f} 0.000 {charge:.6f}\n"

    # Write CIF content with charges to a file
    output_cif_with_charges_path = f"{sys}/{sys}_charges.cif"
    with open(output_cif_with_charges_path, 'w') as cif_file:
        cif_file.write(cif_content_with_charges)

    print(f"CIF file generated: {output_cif_with_charges_path}")


for sys in glob('*-MOF-*'):
    xyz_file = open(os.path.join(sys,'DDEC6_even_tempered_net_atomic_charges.xyz')).readlines()
    xyz2cif(xyz_file,sys)
