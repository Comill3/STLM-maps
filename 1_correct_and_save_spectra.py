# -*- coding: utf-8 -*-
"""
Removing the cosmic rays from the spectra and saving the corrected spectra.

Created on Wed Aug 18 15:46:36 2021

@author: Mylène Sauty
"""

import os
import numpy as np
from scipy.interpolate import interp1d
from utils_EA_0 import LoadWavelet1D, CleanLine

DATA_PATH = 'Demo_Data'
DEDICATED_FOLDER = 'STL_8'
DATA_PATH_2 = os.path.join(DATA_PATH, DEDICATED_FOLDER)
DATA_NAME = "280825-8-STL.txt"


def find_cosmic_rays_wavelets(data_line):
    """Find the cosmic rays in the data_line using wavelets.

    Args:
        data_line (np.array): spectral data

    Returns:
        PosRC (list): Positions of the points touched by the cosmic rays
        _ (None): Placeholder for the second return value
        error (float): Error
    """
    # parametres
    eps = 0.03  # 0.01 #0.02 #0.03
    seuil_bas = 575  # 620 #700 #500
    j = 0
    N = 2048  # 256 #512 #1024 #2048

    # Loading the wavelet
    # Generate the fourier transform of a bump steerable wavelet at a 2^j
    # characteristic scale, defined on an array of size N
    wave1 = LoadWavelet1D(j, N)

    # Find points touched by the cosmic rays
    PosRC, _, error = CleanLine(data_line, wave1, eps, seuil_bas)

    return PosRC, error


def find_points_to_remove(PosRC):
    """
        Adds an index after each cosmic ray position to PosRC.

        Args:
            PosRC (list): Positions of the points touched by the cosmic rays

        Returns:
            indices (list): Updated list of positions including the index after each cosmic ray
    """

    indices = []

    for i in PosRC:
        indices.append(i)
        indices.append(i + 1)
    indices = list(set(indices))
    indices.sort()

    return indices


def remove_cosmic_rays_2(intensity, wavelength, PosRC, BGfix):
    """ Removes cosmic rays from the intensity data by interpolating over the affected points.

    Args:
        intensity (np.array): Light intensity array.
        wavelength (np.array): Wavelength array corresponding to the intensity values.
        PosRC (list): List of positions (indices) of the points affected by cosmic rays.
        BGfix (float): Background intensity value to use for correction at the boundaries.

    Returns:
        tuple: A tuple containing:
            - intensity_corrected (np.array): The corrected intensity 
            array with cosmic rays removed.
            - n (int): The number of points that were corrected.
    """
    # PosRC, error = find_cosmic_rays_wavelets(Intensity)
    indices = PosRC

    # print(indices)
    n = len(indices)
    k = len(intensity)

    intensity_corrected = intensity.copy()

    if 0 in indices:
        indices.remove(0)
        intensity_corrected[0] = BGfix
    if k - 1 in indices:
        indices.remove(k - 1)
        intensity_corrected[k - 1] = BGfix
    if k in indices:
        indices.remove(k)

    wavelength_2 = np.delete(wavelength, indices)
    intensity2 = np.delete(intensity_corrected, indices)

    f2 = interp1d(wavelength_2, intensity2, kind='linear')

    for i in indices:
        intensity_corrected[i] = f2(wavelength[i])

    return intensity_corrected, n


def bin_array(array):
    """ Reduces the size of the input array by averaging every four elements.

    Args:
        array (numpy.ndarray): The input array to be binned. The length of the array should be a multiple of 4.

    Returns:
        numpy.ndarray: A new array where each element is the average of four elements from the input array.
    """

    l = np.shape(array)[0] // 4
    new_array = np.zeros(l, dtype=array.dtype)

    for i in range(l):
        new_array[i] = (array[4 * i] + array[4 * i + 1] +
                        array[4 * i + 2] + array[4 * i + 3]) / 4

    return new_array


def read_and_correct_data(datapath, dataname):
    """Reads the data from the file, corrects the cosmic rays and saves the corrected spectra.

    Args:
        datapath (str): data path
        dataname (str): data name
    """
    BGfix = 600
    Binning = False

    # Load data
    STL_data = np.loadtxt(os.path.join(datapath, dataname), skiprows=11).T
    (m, l) = np.shape(STL_data)
    print(m, l)

    # Get data from file
    wavelength = STL_data[1][:]  # [:-4]
    energy = 1240 / wavelength

    if Binning:
        wavelength = bin_array(wavelength)
        energy = 1240 / wavelength

    spectra = STL_data[2:m]

    name = 'Wavelength'
    if Binning:
        name += '_Binning'

    wavelength = np.array(wavelength)
    np.save(os.path.join(datapath, name), wavelength)

    # Select only spectra without cosmic rays to apply the pca
    corrected_spectra = []
    selected_spectra = []

    for j in range(0, m - 2, 1):
        Intj = spectra[j][:]  # [:-4]

        if Binning:
            Intj = bin_array(Intj)

        # Remove cosmic ray from data (interpolate with neighbours to replace the point) and remove constant background
        PosRC, error = find_cosmic_rays_wavelets(Intj)
        PosRC = find_points_to_remove(PosRC)
        Intensity_corrected, n = remove_cosmic_rays_2(
            Intj, wavelength, PosRC, BGfix)
        if n != 0:
            PosRC2, error = find_cosmic_rays_wavelets(Intensity_corrected)
            PosRC2 = find_points_to_remove(PosRC2)
            PosRC = list(set(PosRC + PosRC2))
            Intensity_corrected, n1 = remove_cosmic_rays_2(
                Intensity_corrected, wavelength, PosRC, BGfix)

        n = len(PosRC)

        Intj_corrected = Intensity_corrected - BGfix

        corrected_spectra.append(Intj_corrected)
        if n == 0:
            selected_spectra.append(Intj_corrected)

        if j % 10 == 0:
            print(j)

    name = 'Corrected_spectra'
    if Binning == True:
        name += '_Binning'

    name2 = 'Selected_spectra_for_PCA'
    if Binning == True:
        name2 += '_Binning'

    corrected_spectra = np.array(corrected_spectra)
    np.save(os.path.join(datapath, name), corrected_spectra)
    selected_spectra = np.array(selected_spectra)
    np.save(os.path.join(datapath, name2), selected_spectra)

read_and_correct_data(DATA_PATH_2, DATA_NAME)
