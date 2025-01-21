import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from Fit_functions import get_data_for_fit, gaussian

datapath = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-03 #SG19 ALE W UHV LN/STL_37'
PCA = False
smooth = False
n = 7
Eph = 0.09  # Phonon replica energy eV
BsA = 0.2


# liste_high = [31, 95, 156, 186, 217, 218, 249, 375, 517, 518, 519, 520, 548, 549, 550, 551, 552, 553, 555, 581, 582, 583, 584, 585, 723]
# liste_room = [921, 922, 923, 924]

# liste_all = list(range(0, 1024))
# liste_exclu = liste_high
# liste_other = [x for x in liste_all if x not in liste_exclu]

STL_data, STL_data_init, Wavelength, Energy, savepath, m, l = get_data_for_fit(PCA, datapath, n)

foldcontents = os.listdir(savepath)
foldername = 'Fit_gaussian_plot'
if foldername not in foldcontents:
    os.makedirs(os.path.join(savepath, foldername))
savingpath = os.path.join(savepath, foldername)

couleur = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
           'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'tab:magenta']

start = 0
stop = m
coeff_fit_8G = []

W1 = 1240 / 1.5
W2 = 1240 / 3.5
i1 = np.argmin(np.abs(Wavelength - W1))
i2 = np.argmin(np.abs(Wavelength - W2))
#    i1 = 0
#    i2 = -1

liste_spectra = np.arange(0, m, 1)

for j in liste_spectra:
            
        #if j in liste_other:

            def func1(x, A1, B1, C1, E1, sigma1):
                return gaussian(x, np.abs(A1), E1, sigma1) + gaussian(x, np.abs(B1), E1 - Eph, sigma1) + gaussian(x, np.abs(C1), E1 - 2*Eph, sigma1)

            print(j)

            Spectrum = STL_data[j]
            Spectrum_init = STL_data_init[j]

            if smooth == True:
                Spectrum = savgol_filter(Spectrum, 5, 1)

            ##fit gaussian

            A0 = max(Spectrum) - 2 * np.std(Spectrum[-100:-1])
            E0 = Energy[np.argmax(Spectrum)]

            try:
                popt, pcov = curve_fit(func1, Energy[i1:i2], Spectrum[i1:i2], p0=(A0, A0*BsA, A0*BsA**2, 2.3, 0.05), maxfev=100000)
                A1, B1, C1, E1, sigma1  = popt

                coeff_fit_8G.append([8, np.abs(A1), E1, np.abs(sigma1), np.abs(B1), E1 - Eph, np.abs(sigma1), 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

                fit = func1(Energy, A1, B1, C1, E1, sigma1)


            except RuntimeError:
                  
                coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

            
            if A1 < 1:
                coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

            # fig = plt.figure(figsize=[7, 5])
            # ax = fig.add_subplot(111)

            # s = f'A1 = {round(A1, 0)}, E1 = {round(E1, 3)}, sigma1 = {round(sigma1, 4)}'
            # s2 = f'B1/A1 = {round(B1/A1, 3)}'
            # plt.plot(Energy[i1:i2], Spectrum_init[i1:i2], 'c', label='raw spectra')
            # if PCA:
            #     plt.plot(Energy[i1:i2], Spectrum[i1:i2], 'gray', label='PCA data')
            # elif smooth:
            #     plt.plot(Energy[i1:i2], Spectrum[i1:i2], 'gray', label='smoothed data')
            # plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(A1), E1, sigma1), color=couleur[1])
            # plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(B1), E1 - Eph, sigma1), color=couleur[2])
            # plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(C1), E1 - 2*Eph, sigma1), color=couleur[3])
            # plt.plot(Energy[i1:i2], fit[i1:i2], 'k', label='fit')

            # plt.text(0.5, 1.07, s, ha='center', va='center', transform=ax.transAxes, fontsize=10)
            # plt.text(0.14, 0.92, s2, ha='center', va='center', transform=ax.transAxes, fontsize=12)
            # plt.rcParams.update({'font.size': 12})
            # plt.legend(loc='upper right')
            # plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
            # ax.tick_params(axis='x', labelsize=14)
            # ax.tick_params(axis='y', labelsize=14)
            # plt.xlabel('Energy (eV)', fontsize=16)
            # plt.ylabel('Intensity (arb. units)', fontsize=16)

            # savename = os.path.join(savingpath, 'Spectrum' + str(j) + '.png')
            # plt.savefig(savename, dpi=300, bbox_inches='tight')

            # plt.show()

        # if j in liste_high:

        #     def func1(x, A1, A2, B1, C1, E1, E2, sigma1, sigma2):
        #         return gaussian(x, np.abs(A1), E1, sigma1) + gaussian(x, np.abs(A2), E2, sigma2) + gaussian(x, np.abs(B1), E1 - Eph, sigma1) + gaussian(x, np.abs(C1), E1 - 2*Eph, sigma1)

        #     print(j)

        #     Spectrum = STL_data[j]
        #     Spectrum_init = STL_data_init[j]

        #     if smooth == True:
        #         Spectrum = savgol_filter(Spectrum, 5, 1)

        #     ##fit gaussian

        #     A0 = max(Spectrum) - 2 * np.std(Spectrum[-100:-1])
        #     E0 = Energy[np.argmax(Spectrum)]

        #     fig = plt.figure(figsize=[7, 5])
        #     ax = fig.add_subplot(111)

        #     try:
        #         popt, pcov = curve_fit(func1, Energy[i1:i2], Spectrum[i1:i2], p0=(A0, A0*0.3, A0*BsA, A0*BsA**2, 2.85, 3, 0.01, 0.01), maxfev=100000)
        #         A1, A2, B1, C1, E1, E2, sigma1, sigma2 = popt

        #         coeff_fit_8G.append([8, np.abs(A1), E1, sigma1, np.abs(B1), E1 - Eph, sigma1, np.abs(A2), E2, sigma2, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

        #         fit = func1(Energy, A1, A2, B1, C1, E1, E2, sigma1, sigma2)

        #     except RuntimeError:
                  
        #         coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, 0, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

        #     s = f'A1 = {round(A1, 0)}, E1 = {round(E1, 3)}, sigma1 = {round(sigma1, 4)}'
        #     #s2 = 'B1/A1 = {}'.format(round(B1/A1, 3))
        #     plt.plot(Energy[i1:i2], Spectrum_init[i1:i2], 'c', label='raw spectra')
        #     if PCA:
        #         plt.plot(Energy[i1:i2], Spectrum[i1:i2], 'gray', label='PCA data')
        #     elif smooth:
        #         plt.plot(Energy[i1:i2], Spectrum[i1:i2], 'gray', label='smoothed data')
        #     plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(A1), E1, sigma1), color=couleur[1])
        #     plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(B1), E1 - Eph, sigma1), color=couleur[2])
        #     plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(C1), E1 - 2*Eph, sigma1), color=couleur[3])
        #     plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(A2), E2, sigma2), color=couleur[4])
        #     plt.plot(Energy[i1:i2], fit[i1:i2], 'k', label='fit')

        #     plt.text(0.5, 1.07, s, ha='center', va='center', transform=ax.transAxes, fontsize=10)
        #     #plt.text(0.14, 0.92, s2, ha='center', va='center', transform=ax.transAxes, fontsize=12)
        #     plt.rcParams.update({'font.size': 12})
        #     plt.legend(loc='upper right')
        #     plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
        #     ax.tick_params(axis='x', labelsize=14)
        #     ax.tick_params(axis='y', labelsize=14)
        #     plt.xlabel('Energy (eV)', fontsize=16)
        #     plt.ylabel('Intensity (arb. units)', fontsize=16)

        #     savename = os.path.join(savingpath, 'Spectrum' + str(j) + '.png')
        #     plt.savefig(savename, dpi=300, bbox_inches='tight')

        # #     #plt.show()

        # # if j in liste_high2:

        # #     def func1(x, A1, A2, A3, B1, C1, E1, E2, E3, sigma1, sigma2, sigma3):
        # #         return gaussian(x, np.abs(A1), E1, sigma1) + gaussian(x, np.abs(A2), E2, sigma2) + gaussian(x, np.abs(B1), E1 - Eph, sigma1) + gaussian(x, np.abs(C1), E1 - 2*Eph, sigma1) + gaussian(x, np.abs(A3), E3, sigma3)

        # #     print(j)

        # #     Spectrum = STL_data[j]
        # #     Spectrum_init = STL_data_init[j]

        # #     if smooth == True:
        # #         Spectrum = savgol_filter(Spectrum, 5, 1)

        # #     ##fit gaussian

        # #     A0 = max(Spectrum) - 2 * np.std(Spectrum[-100:-1])
        # #     E0 = Energy[np.argmax(Spectrum)]

        # #     fig = plt.figure(figsize=[7, 5])
        # #     ax = fig.add_subplot(111)

        # #     try:
        # #         popt, pcov = curve_fit(func1, Energy[i1:i2], Spectrum[i1:i2], p0=(A0*0.5, A0, A0, A0*BsA, A0*BsA**2, 2.3, 2.5, 2.6, 0.03, 0.01, 0.01), maxfev=100000)
        # #         A1, A2, A3, B1, C1, E1, E2, E3, sigma1, sigma2, sigma3 = popt

        # #         coeff_fit_8G.append([8, np.abs(A1), E1, sigma1, np.abs(B1), E1 - Eph, sigma1, np.abs(A2), E2, sigma2, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

        # #         fit = func1(Energy, A1, A2, A3, B1, C1, E1, E2, E3, sigma1, sigma2, sigma3)

        # #     except RuntimeError:
                  
        # #         coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, 0, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])

        # #     s = f'A1 = {round(A1, 0)}, E1 = {round(E1, 3)}, sigma1 = {round(sigma1, 4)}'
        # #     #s2 = 'B1/A1 = {}'.format(round(B1/A1, 3))
        # #     plt.plot(Energy[i1:i2], Spectrum_init[i1:i2], 'c', label='raw spectra')
        # #     if PCA:
        # #         plt.plot(Energy[i1:i2], Spectrum[i1:i2], 'gray', label='PCA data')
        # #     elif smooth:
        # #         plt.plot(Energy[i1:i2], Spectrum[i1:i2], 'gray', label='smoothed data')
        # #     plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(A1), E1, sigma1), color=couleur[1])
        # #     plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(B1), E1 - Eph, sigma1), color=couleur[2])
        # #     plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(A2), E2, sigma2), color=couleur[4])
        # #     plt.plot(Energy[i1:i2], gaussian(Energy[i1:i2], np.abs(A3), E3, sigma3), color=couleur[5])
        # #     plt.plot(Energy[i1:i2], fit[i1:i2], 'k', label='fit')

        # #     plt.text(0.5, 1.07, s, ha='center', va='center', transform=ax.transAxes, fontsize=10)
        # #     #plt.text(0.14, 0.92, s2, ha='center', va='center', transform=ax.transAxes, fontsize=12)
        # #     plt.rcParams.update({'font.size': 12})
        # #     plt.legend(loc='upper right')
        # #     plt.legend(bbox_to_anchor=(1, 1), fontsize=14)
        # #     ax.tick_params(axis='x', labelsize=14)
        # #     ax.tick_params(axis='y', labelsize=14)
        # #     plt.xlabel('Energy (eV)', fontsize=16)
        # #     plt.ylabel('Intensity (arb. units)', fontsize=16)

        # #     savename = os.path.join(savingpath, 'Spectrum' + str(j) + '.png')
        # #     plt.savefig(savename, dpi=300, bbox_inches='tight')

        # # #     #plt.show()

        # if j in liste_room:

        #     coeff_fit_8G.append([8, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, 0, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan, 0, np.nan, np.nan])
        


coeff_fit_8G = np.array(coeff_fit_8G)

np.savetxt(os.path.join(savepath, 'coeff_8Gaussian_fit_from{}to{}.txt'.format(start, stop)), coeff_fit_8G, comments='# number of gaussian, and all gaussian param in order ampl, Energy, sigma', delimiter='\t')