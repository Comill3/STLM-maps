import os
import numpy as np
from Images_functions import create_grid, get_concat_h, get_concat_v, join_images

grid_size = [32,32]
image_size = 500 # nm
step_linescan = image_size/grid_size[1]
missing_spectra = []
perverted_spectra = [i+1 for i in missing_spectra]
missing_spectra = np.concatenate([missing_spectra, perverted_spectra])
data_path = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-03 #SG19 ALE W UHV LN/STL_40/Plot_spectra'
grid_path = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-03 #SG19 ALE W UHV LN/STL_40/Plot_spectra/grid'



foldcontents = os.listdir(data_path)
if 'grid' not in foldcontents:
    os.makedirs(os.path.join(data_path,'grid'))
saving_path = os.path.join(data_path,'grid')

# # Creation of the grids
# create_grid(grid_size, image_size, step_linescan, missing_spectra, saving_path)

# Concatenate spectra and grids
name = '060924-40-STLvsEn_waterfall'
# name = '060924-40-STLvsEn_2D'
number_line = grid_size[0]
mode = 'horizontal' #'vertical' #'horizontal'

join_images(data_path, grid_path, mode, number_line, name, grid_size)
