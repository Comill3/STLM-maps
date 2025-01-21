import os
import numpy as np

def gaussian(x, A, E, sigma):
    """Gaussian function

    Args:
        x (array): x data
        A (float): amplitude
        E (float): energy position
        sigma (float): 2,355*sigma = FWHM

    Returns:
        array: A * np.exp(-(x - E) ** 2 / (2 * sigma ** 2))
    """
    return A * np.exp(-(x - E) ** 2 / (2 * sigma ** 2))

def get_data_for_fit(PCA, datapath, n):
    """_summary_

    Args:
        PCA (boolean): get the PCA or raw data
        datapath (string): data path
        n (int): number of eigen vectors

    Returns:
        _type_: _description_
    """

    # Get the PCA data
    if PCA == True:
        name_folder = '{}vecteurs_propres'.format(n)
        basispath = os.path.join(datapath, name_folder)
        basisname = "{}vecteurs_propres.npy".format(n)
        dataname = 'coeff_projection_from0to1024.npy'

        # Get data projected on PCA basis
        coeff = np.load(os.path.join(basispath, dataname))
        vectors = np.load(os.path.join(basispath, basisname))
        for i in range(np.shape(vectors)[0]):
            vectors[i] = vectors[i] / max(vectors[i])
        # Data projected on PCA basis
        STL_data = np.matmul(coeff, vectors)

        (m, l) = np.shape(STL_data)
        print(m, l)

        # Get already corrected data
        name_data = 'Corrected_spectra'
        name_wl = 'Wavelength'
        # Corrected data
        STL_data_init = np.load(os.path.join(datapath, name_data + '.npy'))
        Wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        Energy = 1240 / Wavelength
        savepath = basispath

    # Or the raw ones
    else:
        # Get already corrected data
        name_data = 'Corrected_spectra'  # 'Interference_Corrected_spectra' #'Corrected_spectra'
        name_wl = 'Wavelength'
        STL_data = np.load(os.path.join(datapath, name_data + '.npy'))
        Wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        Energy = 1240 / Wavelength

        (m, l) = np.shape(STL_data)
        print(m, l)

        STL_data_init = np.load(os.path.join(datapath, name_data + '.npy'))
        Wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
        Energy = 1240 / Wavelength

        savepath = datapath

    return STL_data, STL_data_init, Wavelength, Energy, savepath, m, l

def get_gaussian_coeff(PCA, datapath, n, Gnumber, start, stop):
    """Extract gaussian coeff from the text file

    Args:
        PCA (boolean): PCA data ?
        datapath (string): data path
        n (int): number of eigen vectors
        Gnumber (string): Number of gaussian in the array
        start (int): 0
        stop (int): number of spectra - 1

    Returns:
        string: name of the file and the saving path
    """

    # PCA data
    if PCA == True:
        name_folder = '{}vecteurs_propres'.format(n)
        basispath = os.path.join(datapath, name_folder)
        dataname = 'coeff_' + Gnumber + 'Gaussian_fit_from{}to{}.txt'.format(start, stop)

    else:
        basispath = datapath
        dataname = 'coeff_' + Gnumber + 'Gaussian_fit_from{}to{}.txt'.format(start, stop)

    return dataname, basispath