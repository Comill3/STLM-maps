"""Plot the STL spectra in a grid and the sum of the spectra in the grid, using PCA basis vectors."""

import os

import numpy as np

from plot_functions import plot_stl_in_grid_PCA, plot_sum_spectra_PCA

datapath = "Demo_Data"

n = 7  # Number of eigen vectors
m = 1024  # Number of spectra

# Plot STL spectra in a grid, with one file in one dedicated folder
dedicated_folder = "STL_8"
datapath2 = os.path.join(datapath, dedicated_folder)
stl_to_plot = "250825-8-STL.txt"

Smooth = False

v_min = 1.5
v_max = 3.5
Full = False

name_folder = f"{n}vecteurs_propres"
basispath = os.path.join(datapath2, name_folder)
basisname = f"{n}vecteurs_propres.npy"
dataname = f"coeff_projection_from0to{m}.npy"

path = [datapath2, basispath, basisname, dataname]

step_linescan = 1000 / 32  # in nm #for STL measurements
grid_size = 32  # grid
missing_spectra = []
perverted_spectra = [i + 1 for i in missing_spectra]

# Dictionnary
param_dict = {}
data_param = np.loadtxt(os.path.join(datapath, "Data_parameters.txt"))
(p, q) = np.shape(data_param)
param_dict = {}
for i in range(p):
    param_dict[data_param[i, 0]] = (
        data_param[i, 1],
        data_param[i, 2],
        data_param[i, 3],
    )

plot_stl_in_grid_PCA(
    stl_to_plot,
    path,
    param_dict,
    step_linescan,
    grid_size,
    missing_spectra,
    perverted_spectra,
    Full,
    v_min,
    v_max,
    Smooth,
    SVG=False,
)

plot_sum_spectra_PCA(
    stl_to_plot,
    path,
    param_dict,
    missing_spectra,
    perverted_spectra,
    Full,
    v_min,
    v_max,
    Smooth,
)
