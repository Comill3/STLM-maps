"""Find the local maxima of each spectra."""

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter
from fit_functions import get_data_for_fit

DATA_PATH = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2023-06-19 #SG19 ref(1) W UHV RT/STL_35'
PCA = True
n = 13

stl_data, stl_data_init, wavelength, energy, savepath, m, l = get_data_for_fit(PCA, DATA_PATH, n)

foldcontents = os.listdir(savepath)
foldername = 'Local_maxima_plot'
if foldername not in foldcontents:
    os.makedirs(os.path.join(savepath, foldername))
savingpath = os.path.join(savepath, foldername)

couleur = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
           'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'tab:magenta']

start = 0
stop = 1024
coeff_fit_8G = []

W1 = 1240 / 1.5
W2 = 1240 / 3.5
i1 = np.argmin(np.abs(wavelength - W1))
i2 = np.argmin(np.abs(wavelength - W2))
#    i1 = 0
#    i2 = -1

liste_spectra = np.arange(0, 1024, 1)

# Changing parameters
dE_0 = (energy[-1] - energy[0])/len(energy) #eV per point
dE = 0.05 # min distance between peaks (eV)
distance_min = int(dE/dE_0)  # min distance between peaks (points)
prominence_min = 0.01  # relative prominence (0-1)

liste_totale = []

for j in liste_spectra:

    liste_peaks = []

    print(j)
        
    spectrum = stl_data[j]
    spectrum_smooth = savgol_filter(spectrum, 20, 1)
    spectrum_init = stl_data_init[j]

    peaks, prop = find_peaks(spectrum_smooth[i1:i2], 
                             distance=distance_min, 
                             prominence=prominence_min*np.max(spectrum[i1:i2]),
                             height=20)
    
    s = f'# {len(peaks)}'

    fig, ax = plt.subplots()
    plt.text(0.14, 0.92, s, ha='center', va='center', transform=ax.transAxes, fontsize=12)
    plt.plot(energy[i1:i2], spectrum_init[i1:i2], color='c', label='Raw data')
    plt.plot(energy[i1:i2], spectrum[i1:i2], color='grey', label='Denoised data')
    plt.plot(energy[i1:i2], spectrum_smooth[i1:i2], color='k', label='Denoised and smoothed data')
    plt.plot(energy[peaks], spectrum_smooth[peaks], 'rx', label='Identifided peaks')
    plt.legend()
    plt.xlabel('Energy (eV)', fontsize=12)
    plt.ylabel('Intensity (arb. units)', fontsize=12)
    # plt.show()

    for i in energy[peaks]:
        if i<2 or spectrum_smooth[np.argmin(np.abs(energy - i))]<20:
            liste_peaks.append((np.nan, np.nan))
        else:
            liste_peaks.append((i, spectrum_smooth[np.argmin(np.abs(energy - i))]))

    num_valid_peaks = sum(1 for peak in liste_peaks if not np.isnan(peak[0]))
    liste_totale.append((num_valid_peaks, liste_peaks))
    
    
    save_name = os.path.join(savingpath, f'Spectrum_{j}.png')
    plt.savefig(save_name, dpi=300, bbox_inches='tight')

with open(os.path.join(savepath, 'local_maxima.pkl'), 'wb') as f:
    pickle.dump(liste_totale, f)

# np.save(os.path.join(savepath, f'list_peaks{start}to{stop}.txt'), liste_peaks, comments='# Local maxima position and amplitude', delimiter='\t')