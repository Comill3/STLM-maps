"""Create grid. Concatenate the spectra and the grids."""

import os
import numpy as np
from images_functions import create_grid, join_images

grid_size = [32,32]
image_size = 1000 # nm
step_linescan = image_size/grid_size[1]
missing_spectra = [55, 260, 412, 429, 612]
perverted_spectra = [i+1 for i in missing_spectra]
missing_spectra = np.concatenate([missing_spectra, perverted_spectra])
data_path = 'Demo_Data/STL_8/7vecteurs_propres/Plot_spectra'
grid_path = 'Demo_Data/STL_8/7vecteurs_propres/Plot_spectra/grid'


foldcontents = os.listdir(data_path)
if 'grid' not in foldcontents:
    os.makedirs(os.path.join(data_path,'grid'))
saving_path = os.path.join(data_path,'grid')

# Creation of the grids
create_grid(grid_size, image_size, step_linescan, missing_spectra, saving_path, SVG=False)

### Uncomment when the grids are made ###

# # Concatenate spectra and grids
# name = '190623-35-STLvsEn_waterfall'
# # name = '190623-35-STLvsEn_2D'
# number_line = grid_size[0]
# mode = 'horizontal' #'vertical' #'horizontal'

# join_images(data_path, grid_path, mode, number_line, name, grid_size)
