"""Specific function used to fit STLM data."""

import os
import numpy as np

def gaussian(x, A, E, sigma):
    """Gaussian function.

    Args:
        x (array): x data.
        A (float): Amplitude.
        E (float): Energy position.
        sigma (float): 2,355*sigma = FWHM

    Returns:
        array: Gaussian function.
    """
    return A * np.exp(-(x - E) ** 2 / (2 * sigma ** 2))

def get_data_for_fit(PCA, datapath, n):
    """Get the spectral data for the fit.

    Args:
        PCA (bool): If True get the PCA data, else the raw/smoothed data.
        datapath (str): Data path.
        n (int): Number of eigen vectors.

    Returns:
        Array: STL data.
        Array: initial STL data.
        Array: Wavelength.
        Array: Energy.
        String: Save path.
        Int: Number of X pixels (1600).
        Int: Number of colomns, number of spectra.
    """

    # Get the PCA data
    if PCA:
        name_folder = f'{n}vecteurs_propres'
        basispath = os.path.join(datapath, name_folder)
        basisname = f'{n}vecteurs_propres.npy'
        dataname = 'coeff_projection_from0to1024.npy'

        # Get data projected on PCA basis
        coeff = np.load(os.path.join(basispath, dataname))
        vectors = np.load(os.path.join(basispath, basisname))
        for i in range(np.shape(vectors)[0]):
            vectors[i] = vectors[i] / max(vectors[i])
        # Data projected on PCA basis
        stl_data = np.matmul(coeff, vectors)

        (m, l) = np.shape(stl_data)
        print(m, l)

        # Get already corrected data
        name_data = 'Corrected_spectra'
        name_wl = 'Wavelength'
        # Corrected data
        stl_data_init = np.load(os.path.join(datapath, name_data + '.npy'))
        wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        energy = 1240 / wavelength
        savepath = basispath

    # Or the raw ones
    else:
        # Get already corrected data
        name_data = 'Corrected_spectra'  # 'Interference_Corrected_spectra' #'Corrected_spectra'
        name_wl = 'Wavelength'
        stl_data = np.load(os.path.join(datapath, name_data + '.npy'))
        wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        energy = 1240 / wavelength

        (m, l) = np.shape(stl_data)
        print(m, l)

        stl_data_init = np.load(os.path.join(datapath, name_data + '.npy'))
        wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        energy = 1240 / wavelength

        savepath = datapath

    return stl_data, stl_data_init, wavelength, energy, savepath, m, l

def get_gaussian_coeff(PCA, datapath, n, G, start, stop):
    """Extract gaussian coeff from the text file to plot the maps.

    Args:
        PCA (bool): If True get the PCA data, else the raw/smoothed data.
        datapath (str): Data path.
        n (int): Number of eigen vectors.
        G (str): Number of gaussian saved in the array.
        start (int): 0.
        stop (int): Number of spectra - 1.

    Returns:
        String: Name of the file
        String: Basis path.
    """

    # PCA data
    if PCA:
        name_folder = f'{n}vecteurs_propres'
        basispath = os.path.join(datapath, name_folder)
        dataname = 'coeff_' + G + f'Gaussian_fit_from{start}to{stop}.txt'

    else:
        basispath = datapath
        dataname = 'coeff_' + G + f'Gaussian_fit_from{start}to{stop}.txt'

    return dataname, basispath

def get_data_for_max(PCA, datapath, n):
    """Get the spectral data to find the position of the maximum.

    Args:
        PCA (bool): If True get the PCA data, else the raw/smoothed data.
        datapath (str): Data path.
        n (int): Number of eigen vectors.

    Returns:
        Array: STL data.
        Array: initial STL data.
        Array: Wavelength.
        Array: Energy.
        String: Save path.
        Int: Number of X pixels (1600).
        Int: Number of colomns, number of spectra.
    """

    # Get the PCA data
    if PCA:
        name_folder = f'{n}vecteurs_propres'
        basispath = os.path.join(datapath, name_folder)
        basisname = f'{n}vecteurs_propres.npy'
        dataname = 'coeff_projection_from0to1024.npy'

        # Get data projected on PCA basis
        coeff = np.load(os.path.join(basispath, dataname))
        vectors = np.load(os.path.join(basispath, basisname))
        for i in range(np.shape(vectors)[0]):
            vectors[i] = vectors[i] / max(vectors[i])
        # Data projected on PCA basis
        stl_data = np.matmul(coeff, vectors)

        (m, l) = np.shape(stl_data)
        print(m, l)

        # Get already corrected data
        name_data = 'Corrected_spectra'
        name_wl = 'Wavelength'
        # Corrected data
        stl_data_init = np.load(os.path.join(datapath, name_data + '.npy'))
        wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        energy = 1240 / wavelength
        savepath = basispath

    # Or the raw ones
    else:
        # Get already corrected data
        name_data = 'Corrected_spectra'  # 'Interference_Corrected_spectra' #'Corrected_spectra'
        name_wl = 'Wavelength'
        stl_data = np.load(os.path.join(datapath, name_data + '.npy'))
        wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        energy = 1240 / wavelength

        stl_data_init = np.load(os.path.join(datapath, name_data + '.npy'))
        wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        energy = 1240 / wavelength

    return stl_data, stl_data_init, wavelength, energy