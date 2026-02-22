"""Fit the STLM spectra with a combination of gaussians."""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from fit_functions import get_data_for_fit, gaussian

DATA_PATH = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2025-08-25 #SG19 ALE UHV W RT/STL_13'
PCA = True
Smooth = False
n = 5
Eph = 0.09  # Phonon replica energy in eV
BsA = 0.2 #Ration between QW emission and phonon replica

stl_data, stl_data_init, wavelength, energy, savepath, m, l = get_data_for_fit(PCA, DATA_PATH, n)

foldcontents = os.listdir(savepath)
foldername = 'Fit_gaussian_plot'
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

# liste_2_5 = [6, 38, 70, 99, 102, 103, 136, 193, 197, 199, 224, 225, 230, 403, 414, 703]
# liste_2_6 = [622, 623, 624, 625, 626, 654, 655, 656, 657]
# list_both = [192, 223]

# # liste_spe = liste_2_5 + liste_2_6 + list_both
# liste_all = [i for i in liste_spectra if i not in liste_2_6]

for j in liste_spectra:

    # if j in liste_all:
        def func1(x, A1, B1, E1, sigma1):
            """Fitting function using a combination of Gaussian: A1, B1, C1, A2, B2,
            """
            return gaussian(x, np.abs(A1), E1, sigma1) + gaussian(x, np.abs(B1), E1 - Eph, sigma1)

        print(j)

        spectrum = stl_data[j]
        spectrum_init = stl_data_init[j]

        if Smooth:
            spectrum = savgol_filter(spectrum, 5, 1)

        ##fit gaussian

        A0 = max(spectrum) - 2 * np.std(spectrum[-100:-1])
        E0 = energy[np.argmax(spectrum)]

        try:
            popt, pcov = curve_fit(func1, energy[i1:i2], spectrum[i1:i2], p0=(A0, A0*BsA, 2.3, 0.05), maxfev=100000)
            A1, B1, E1, sigma1  = popt

            coeff_fit_8G.append([8, np.abs(A1), E1, np.abs(sigma1), np.abs(B1), E1 - Eph, np.abs(sigma1), 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

            if A1 < 1:
                coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

            fit = func1(energy, A1, B1, E1, sigma1)

        except RuntimeError:  
            coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

        fig = plt.figure(figsize=[7, 5])
        ax = fig.add_subplot(111)

        s = f'A1 = {round(A1, 0)}, E1 = {round(E1, 3)}, sigma1 = {round(sigma1, 4)}'
        s2 = f'B1/A1 = {round(B1/A1, 3)}'
        plt.plot(energy[i1:i2], spectrum_init[i1:i2], 'c', label='raw spectra')
        if PCA:
            plt.plot(energy[i1:i2], spectrum[i1:i2], 'gray', label='PCA data')
        elif Smooth:
            plt.plot(energy[i1:i2], spectrum[i1:i2], 'gray', label='smoothed data')

        plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(A1), E1, sigma1), color=couleur[1])
        plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(B1), E1 - Eph, sigma1), color=couleur[2])
        plt.plot(energy[i1:i2], fit[i1:i2], 'k', label='fit')

        plt.text(0.5, 1.07, s, ha='center', va='center', transform=ax.transAxes, fontsize=10)
        plt.text(0.14, 0.92, s2, ha='center', va='center', transform=ax.transAxes, fontsize=12)
        plt.rcParams.update({'font.size': 12})
        plt.legend(loc='upper right')
        plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
        ax.tick_params(axis='x', labelsize=14)
        ax.tick_params(axis='y', labelsize=14)
        plt.xlabel('Energy (eV)', fontsize=16)
        plt.ylabel('Intensity (arb. units)', fontsize=16)

        savename = os.path.join(savingpath, 'Spectrum' + str(j) + '.png')
        plt.savefig(savename, dpi=300, bbox_inches='tight')

        # plt.show()

    # if j in liste_2_5:

    #     def func1(x, A1, B1, C1, A2, E1, E2, sigma1, sigma2):
    #         """Fitting function using a combination of Gaussian: A1, B1, C1, A2, B2,
    #         """
    #         return gaussian(x, np.abs(A1), E1, sigma1) + gaussian(x, np.abs(B1), E1 - Eph, sigma1) + gaussian(x, np.abs(C1), E1 - 2*Eph, sigma1) + gaussian(x, np.abs(A2), E2, sigma2)

    #     print(j)

    #     spectrum = stl_data[j]
    #     spectrum_init = stl_data_init[j]

    #     if Smooth:
    #         spectrum = savgol_filter(spectrum, 5, 1)

    #     ##fit gaussian

    #     A0 = max(spectrum) - 2 * np.std(spectrum[-100:-1])
    #     E0 = energy[np.argmax(spectrum)]

    #     try:
    #         popt, pcov = curve_fit(func1, energy[i1:i2], spectrum[i1:i2], p0=(A0, A0*BsA, A0*BsA**2, A0*0.5, 2.3, 2.5, 0.05, 0.05), maxfev=100000)
    #         A1, B1, C1, A2, E1, E2, sigma1, sigma2  = popt

    #         coeff_fit_8G.append([8, np.abs(A1), E1, np.abs(sigma1), np.abs(B1), E1 - Eph, np.abs(sigma1), 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #         if A1 < 1:
    #             coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #         fit = func1(energy, A1, B1, C1, A2, E1, E2, sigma1, sigma2)

    #     except RuntimeError:  
    #         coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #     fig = plt.figure(figsize=[7, 5])
    #     ax = fig.add_subplot(111)

    #     s = f'A1 = {round(A1, 0)}, E1 = {round(E1, 3)}, sigma1 = {round(sigma1, 4)}'
    #     s2 = f'B1/A1 = {round(B1/A1, 3)}'
    #     plt.plot(energy[i1:i2], spectrum_init[i1:i2], 'c', label='raw spectra')
    #     if PCA:
    #         plt.plot(energy[i1:i2], spectrum[i1:i2], 'gray', label='PCA data')
    #     elif Smooth:
    #         plt.plot(energy[i1:i2], spectrum[i1:i2], 'gray', label='smoothed data')

    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(A1), E1, sigma1), color=couleur[1])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(B1), E1 - Eph, sigma1), color=couleur[2])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(C1), E1 - 2*Eph, sigma1), color=couleur[3])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(A2), E2, sigma2), color=couleur[4])
    #     plt.plot(energy[i1:i2], fit[i1:i2], 'k', label='fit')

    #     plt.text(0.5, 1.07, s, ha='center', va='center', transform=ax.transAxes, fontsize=10)
    #     plt.text(0.14, 0.92, s2, ha='center', va='center', transform=ax.transAxes, fontsize=12)
    #     plt.rcParams.update({'font.size': 12})
    #     plt.legend(loc='upper right')
    #     plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
    #     ax.tick_params(axis='x', labelsize=14)
    #     ax.tick_params(axis='y', labelsize=14)
    #     plt.xlabel('Energy (eV)', fontsize=16)
    #     plt.ylabel('Intensity (arb. units)', fontsize=16)

    #     savename = os.path.join(savingpath, 'Spectrum' + str(j) + '.png')
    #     plt.savefig(savename, dpi=300, bbox_inches='tight')

    #     # plt.show()

    # if j in liste_2_6:

    #     def func1(x, A1, B1, A2, E1, E2, sigma1, sigma2):
    #         return gaussian(x, np.abs(A1), E1, sigma1) + gaussian(x, np.abs(B1), E1 - Eph, sigma1) + gaussian(x, np.abs(A2), E2, sigma2)

    #     print(j)

    #     spectrum = stl_data[j]
    #     spectrum_init = stl_data_init[j]

    #     if Smooth:
    #         spectrum = savgol_filter(spectrum, 5, 1)

    #     ##fit gaussian

    #     A0 = max(spectrum) - 2 * np.std(spectrum[-100:-1])
    #     E0 = energy[np.argmax(spectrum)]

    #     try:
    #         popt, pcov = curve_fit(func1, energy[i1:i2], spectrum[i1:i2], p0=(A0, A0*BsA, A0*0.5, 2.2, 2.7, 0.05, 0.05), maxfev=100000)
    #         A1, B1, A2, E1, E2, sigma1, sigma2  = popt

    #         coeff_fit_8G.append([8, np.abs(A1), E1, np.abs(sigma1), np.abs(B1), E1 - Eph, np.abs(sigma1), 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #         if A1 < 1:
    #             coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #         fit = func1(energy, A1, B1, A2, E1, E2, sigma1, sigma2)

    #     except RuntimeError:  
    #         coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #     fig = plt.figure(figsize=[7, 5])
    #     ax = fig.add_subplot(111)

    #     s = f'A1 = {round(A1, 0)}, E1 = {round(E1, 3)}, sigma1 = {round(sigma1, 4)}'
    #     s2 = f'B1/A1 = {round(B1/A1, 3)}'
    #     plt.plot(energy[i1:i2], spectrum_init[i1:i2], 'c', label='raw spectra')
    #     if PCA:
    #         plt.plot(energy[i1:i2], spectrum[i1:i2], 'gray', label='PCA data')
    #     elif Smooth:
    #         plt.plot(energy[i1:i2], spectrum[i1:i2], 'gray', label='smoothed data')

    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(A1), E1, sigma1), color=couleur[1])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(B1), E1 - Eph, sigma1), color=couleur[2])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(A2), E2, sigma2), color=couleur[5])
    #     plt.plot(energy[i1:i2], fit[i1:i2], 'k', label='fit')

    #     plt.text(0.5, 1.07, s, ha='center', va='center', transform=ax.transAxes, fontsize=10)
    #     plt.text(0.14, 0.92, s2, ha='center', va='center', transform=ax.transAxes, fontsize=12)
    #     plt.rcParams.update({'font.size': 12})
    #     plt.legend(loc='upper right')
    #     plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
    #     ax.tick_params(axis='x', labelsize=14)
    #     ax.tick_params(axis='y', labelsize=14)
    #     plt.xlabel('Energy (eV)', fontsize=16)
    #     plt.ylabel('Intensity (arb. units)', fontsize=16)

    #     savename = os.path.join(savingpath, 'Spectrum' + str(j) + '.png')
    #     plt.savefig(savename, dpi=300, bbox_inches='tight')

    #     # plt.show()

    # if j in list_both:

    #     def func1(x, A1, B1, C1, A2, A3, E1, E2, E3, sigma1, sigma2, sigma3):
    #         """Fitting function using a combination of Gaussian: A1, B1, C1, A2, B2,
    #         """
    #         return gaussian(x, np.abs(A1), E1, sigma1) + gaussian(x, np.abs(B1), E1 - Eph, sigma1) + gaussian(x, np.abs(C1), E1 - 2*Eph, sigma1) + gaussian(x, np.abs(A2), E2, sigma2) + gaussian(x, np.abs(A3), E3, sigma3)

    #     print(j)

    #     spectrum = stl_data[j]
    #     spectrum_init = stl_data_init[j]

    #     if Smooth:
    #         spectrum = savgol_filter(spectrum, 5, 1)

    #     ##fit gaussian

    #     A0 = max(spectrum) - 2 * np.std(spectrum[-100:-1])
    #     E0 = energy[np.argmax(spectrum)]

    #     try:
    #         popt, pcov = curve_fit(func1, energy[i1:i2], spectrum[i1:i2], p0=(A0, A0*BsA, A0*BsA**2, A0*0.5, A0*0.5, 2.3, 2.5, 2.6, 0.05, 0.05, 0.05), maxfev=100000)
    #         A1, B1, C1, A2, A3, E1, E2, E3, sigma1, sigma2, sigma3  = popt

    #         coeff_fit_8G.append([8, np.abs(A1), E1, np.abs(sigma1), np.abs(B1), E1 - Eph, np.abs(sigma1), 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #         if A1 < 1:
    #             coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #         fit = func1(energy, A1, B1, C1, A2, A3, E1, E2, E3, sigma1, sigma2, sigma3)

    #     except RuntimeError:  
    #         coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

    #     fig = plt.figure(figsize=[7, 5])
    #     ax = fig.add_subplot(111)

    #     s = f'A1 = {round(A1, 0)}, E1 = {round(E1, 3)}, sigma1 = {round(sigma1, 4)}'
    #     s2 = f'B1/A1 = {round(B1/A1, 3)}'
    #     plt.plot(energy[i1:i2], spectrum_init[i1:i2], 'c', label='raw spectra')
    #     if PCA:
    #         plt.plot(energy[i1:i2], spectrum[i1:i2], 'gray', label='PCA data')
    #     elif Smooth:
    #         plt.plot(energy[i1:i2], spectrum[i1:i2], 'gray', label='smoothed data')

    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(A1), E1, sigma1), color=couleur[1])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(B1), E1 - Eph, sigma1), color=couleur[2])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(C1), E1 - 2*Eph, sigma1), color=couleur[3])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(A2), E2, sigma2), color=couleur[4])
    #     plt.plot(energy[i1:i2], gaussian(energy[i1:i2], np.abs(A3), E3, sigma3), color=couleur[5])
    #     plt.plot(energy[i1:i2], fit[i1:i2], 'k', label='fit')

    #     plt.text(0.5, 1.07, s, ha='center', va='center', transform=ax.transAxes, fontsize=10)
    #     plt.text(0.14, 0.92, s2, ha='center', va='center', transform=ax.transAxes, fontsize=12)
    #     plt.rcParams.update({'font.size': 12})
    #     plt.legend(loc='upper right')
    #     plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
    #     ax.tick_params(axis='x', labelsize=14)
    #     ax.tick_params(axis='y', labelsize=14)
    #     plt.xlabel('Energy (eV)', fontsize=16)
    #     plt.ylabel('Intensity (arb. units)', fontsize=16)

    #     savename = os.path.join(savingpath, 'Spectrum' + str(j) + '.png')
    #     plt.savefig(savename, dpi=300, bbox_inches='tight')

    #     # plt.show()

coeff_fit_8G = np.array(coeff_fit_8G)

np.savetxt(os.path.join(savepath, f'coeff_8Gaussian_fit_from{start}to{stop}.txt'), coeff_fit_8G, comments='# number of gaussian, and all gaussian param in order ampl, Energy, sigma', delimiter='\t')
