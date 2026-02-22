""" Plot different maps from the results obtained with '10_local_maxima.py' code. """

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.cm as cm

directory = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2023-06-19 #SG19 ref(1) W UHV RT/STL_35'
PCA = True
n = 13

start = 0
stop = 1024

if PCA:
    folder = f'{n}vecteurs_propres'
    data_path = os.path.join(directory, folder)
else:
    data_path = directory

with open(os.path.join(data_path, 'local_maxima.pkl'), 'rb') as f:
    local_maxima = pickle.load(f)

liste_num = []

for j in np.arange(0, 1024, 1):
    liste_num.append(local_maxima[j][0])

liste_num = np.array(liste_num)

# Plot in grid
image_size = [2000, 2000] # [nm, nm]
grid_size = [32, 32] #32 #grid 
step_linescan = image_size[0]/grid_size[0] #in nm #for STL measurements
missing_spectra = [55, 260, 412, 429, 612]
perverted_spectra = [i+1 for i in missing_spectra]
    
# Add zeros for missing spectra
l2 = len(missing_spectra)
liste_num_modif = liste_num[:stop-start-l2].astype(float)
zeroline = np.nan
for i in missing_spectra:
    if i >= liste_num_modif.shape[0]:
        liste_num_modif = np.append(liste_num_modif, zeroline)
    else:
        liste_num_modif = np.insert(liste_num_modif, i, zeroline)
        
for i in perverted_spectra:
    liste_num_modif[i] = np.nan
    
N = np.reshape(np.copy(liste_num_modif), grid_size)
N = np.flip(np.abs(N), axis = 1)

cmap1 = cm.get_cmap(name='viridis', lut=None)

fig = plt.figure()
plt.imshow(N, cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
clb = plt.colorbar(pad = 0.2)
clb.set_label('# peaks', labelpad=20, y=0.5, rotation=270, fontsize=16)
