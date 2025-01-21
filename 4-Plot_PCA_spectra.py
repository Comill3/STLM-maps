import os
import numpy as np

from plot_functions import plot_STL_in_grid_PCA, plot_sum_spectra_PCA

datapath = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-03 #SG19 ALE W UHV LN'

n = 15  # Number of eigen vectors
m = 1024  # Number of spectra

# Plot STL spectra in a grid, with one file in one dedicated folder
dedicated_folder = 'STL_37'
datapath2 = os.path.join(datapath, dedicated_folder)
STL_to_plot = "060924-37-STL.txt"

smooth = False

Min = 1.5
Max = 3.5
full = False

name_folder = '{}vecteurs_propres'.format(n)
basispath = os.path.join(datapath2, name_folder)
basisname = "{}vecteurs_propres.npy".format(n)
dataname = 'coeff_projection_from0to{}.npy'.format(m)

path = [datapath2, basispath, basisname, dataname]

step_linescan = 1000/32  # in nm #for STL measurements
grid_size = 32  # grid
missing_spectra = []
perverted_spectra = [i+1 for i in missing_spectra] + [332, 333, 334, 335, 336, 337, 338, 339, 340]

param_dict = {}
data_param = np.loadtxt(os.path.join(datapath, "Data_parameters.txt"))
(p, q) = np.shape(data_param)
param_dict = {}
for i in range(p):
    param_dict[data_param[i, 0]] = (
        data_param[i, 1], data_param[i, 2], data_param[i, 3])

plot_STL_in_grid_PCA(STL_to_plot, path, param_dict, step_linescan,
                     grid_size, missing_spectra, perverted_spectra, full, Min, Max, smooth)

plot_sum_spectra_PCA(STL_to_plot, path, param_dict, missing_spectra, perverted_spectra, full, Min, Max, smooth)