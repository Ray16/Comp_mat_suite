import os
import shutil
from glob import glob

DDEC_files_dir = 'DDEC_single_point_files'

for sys in glob('*-MOF-*'):
    for file in glob(DDEC_files_dir+'/*'):
        shutil.copy(file,sys)
