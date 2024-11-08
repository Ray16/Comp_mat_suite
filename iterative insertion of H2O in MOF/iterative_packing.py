import os
import shutil
import subprocess

template_inp = open('template.inp').read()


for insertion_idx in range(100):
    print(f'inserting H2O idx={insertion_idx+1}')
    if insertion_idx == 0:
        current_structure = 'Zr-MOF-801'
    else:
        current_structure = f'h2o_{insertion_idx}'
    num_h2o = 1
    output_f_name = f'h2o_{insertion_idx+1}.inp'
    template_inp_new1 = template_inp.replace('CURRENT_STRUCTURE',current_structure)
    template_inp_new2 = template_inp_new1.replace('OUTPUT_FILE_NAME',output_f_name.split('.inp')[0])
    template_inp_new3 = template_inp_new2.replace('NUM_H2O',str(num_h2o))
    
    with open(output_f_name,'w+') as f:
        f.write(template_inp_new3)
    
    print('Packing h2o...')
    subprocess.run(f'packmol < {output_f_name}',shell=True)
    print('Converting to cif...')
    subprocess.run(f'python xyz2cif.py --idx {insertion_idx+1}',shell=True)
    # optimization
    print('Starting optimization...')
    subprocess.run(f'python MACE_geo_opt_Zr.py --idx {insertion_idx+1}',shell=True)
    # cif2xyz
    subprocess.run(f'python cif2xyz.py --idx {insertion_idx+1}',shell=True)