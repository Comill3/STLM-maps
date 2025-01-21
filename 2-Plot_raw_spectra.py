import os
from plot_functions import *

datapath = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-03 #SG19 ALE W UHV LN'
    
dedicated_folder = 'STL_37'
datapath2 = os.path.join(datapath, dedicated_folder)
STL_to_plot = "050924-30-STL.txt"
Min = 1.5
Max = 3.5
full = False
    
step_linescan = 1000/32 #in nm #for STL measurements
grid_size = 32 #grid 
missing_spectra = [467]
perverted_spectra = [i+1 for i in missing_spectra]
smooth = False

param_dict = {}
data_param = np.loadtxt(os.path.join(datapath, "Data_parameters.txt"))
(p, q) = np.shape(data_param)
param_dict = {}
for i in range(p) : 
    param_dict[data_param[i, 0]] = (data_param[i, 1], data_param[i, 2], data_param[i, 3])

plot_STL_in_grid(STL_to_plot, datapath2, param_dict, step_linescan, grid_size, missing_spectra, perverted_spectra, full, Min, Max, smooth)

plot_sum_spectra(STL_to_plot, datapath2, param_dict, missing_spectra, perverted_spectra, full, Min, Max, smooth)