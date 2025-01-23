"""Plot the raw STLM spectra."""
import os
import numpy as np
from plot_functions import plot_stl_in_grid, plot_sum_spectra

### Data path ###
DATA_PATH = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/' \
            '2024-07-18 #SG19 ALE W UHV LN'
DEDICATED_FOLDER = 'STL_14'
DATA_PATH2 = os.path.join(DATA_PATH, DEDICATED_FOLDER)
STL_to_plot = "190724-14-STL.txt"
v_min = 1.5 #in eV
v_max = 3.5 #in eV
Full = False #True if the full range of the data is to be plotted

### Parameters ###
step_linescan = 2000/32 #in nm #for STL measurements
grid_size = 32 #grid size
missing_spectra = [0, 87, 114, 239, 310, 319, 613, 722, 855, 906]
perverted_spectra = [i+1 for i in missing_spectra]
Smooth = False

### Dictionnary ###
param_dict = {}
data_param = np.loadtxt(os.path.join(DATA_PATH, "Data_parameters.txt"))
(p, q) = np.shape(data_param)
param_dict = {}
for i in range(p): 
    param_dict[data_param[i, 0]] = (data_param[i, 1], data_param[i, 2], data_param[i, 3])

# ### Plotting ###
# plot_STL_in_grid(STL_to_plot, DATA_PATH2, param_dict, step_linescan,
#                  grid_size, missing_spectra, perverted_spectra, Full, v_min, v_max, Smooth)

plot_sum_spectra(STL_to_plot, DATA_PATH2, param_dict, missing_spectra,
                 perverted_spectra, Full, v_min, v_max, Smooth)
