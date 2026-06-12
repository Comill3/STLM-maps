"""Principal Analysis Componants functions for the noise treatement of the acquired data"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

from plot_functions import *

mpl.rcParams["font.family"] = "Arial"


def perform_PCA(datapath, n):
    """Will perfom PCA on the data with the first n eigen vectors

    Args:
        datapath (str): Data path.
        n (int): Number of eigen vectors.
    """

    name_data = "Selected_spectra_for_PCA"  # y data
    name_wl = "Wavelength"  # x data

    selected_spectra = np.load(os.path.join(datapath, name_data + ".npy"))
    print(np.shape(selected_spectra))
    wavelength = np.load(os.path.join(datapath, name_wl + ".npy"))
    energy = 1240 / wavelength

    i1 = 0
    i2 = -1

    # # Specific if spectra (with room light typically) need to be removed
    # list_out = list(range(850, 969))

    # #     if np.max(selected_spectra[i]>10000): # in counts. To be sure this is not from the sample
    # #         list_out.append(i)

    # print(list_out)

    # plt.show()

    # li = []
    # for i in list_out:
    #     plt.plot(selected_spectra[i])
    #     li.append(i)
    # plt.show()

    # print(li)

    # selected_spectra = np.delete(selected_spectra, li, axis=0)

    for i in range(len(selected_spectra)):
        plt.plot(selected_spectra[i])

    plt.show()

    # Perform the PCA on selected spectra
    pca = PCA()
    pca.fit(selected_spectra)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.plot(pca.explained_variance_ratio_[:-1], "kx")
    plt.yscale("log")
    plt.xlim(0, 30)
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)
    # plt.ylim(1e-4, 1)
    plt.xlabel("Vecteur #", fontsize=16)
    plt.ylabel("Rapport de variance", fontsize=16)
    vr = "PCA_variance_ratio"
    plt.savefig(os.path.join(datapath, vr), dpi=300, bbox_inches="tight")
    plt.show()

    # Select only the first n eingenvector
    vecteur_propre = np.zeros((n, np.shape(selected_spectra)[1]))
    for i in range(n):
        vecteur_propre[i] = pca.components_[i]
        # vecteur_propre[i] = gaussian_filter(vecteur_propre[i], sigma=4)
        # vecteur_propre[i] = savgol_filter(vecteur_propre[i], 5, 1)
    plt.plot(energy, vecteur_propre.T)
    plt.xlabel("Energy(eV)")
    plt.ylabel("PCA first eigenvectors")
    plt.grid()

    name_folder = f"{n}vecteurs_propres"

    foldcontents = os.listdir(datapath)
    if name_folder not in foldcontents:
        os.makedirs(os.path.join(datapath, name_folder))
    savingpath = os.path.join(datapath, name_folder)

    plt.savefig(
        os.path.join(savingpath, f"{n}vecteurs_propres.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()

    couleur = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
    ]

    for i in range(n):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(energy[i1:i2], vecteur_propre[i][i1:i2], color=couleur[i])
        plt.xlabel("Energy(eV)", fontsize=16)
        plt.ylabel(f"PCA eigenvector #{i}", fontsize=16)
        ax.tick_params(axis="x", labelsize=16)
        ax.tick_params(axis="y", labelsize=16)
        plt.grid()
        plt.savefig(
            os.path.join(savingpath, f"vecteur_propre_#{i}.png"),
            dpi=300,
            bbox_inches="tight",
        )
        plt.show()

    # Save the orthonormal basis of n first eigenvectors in the data folder
    np.save(os.path.join(savingpath, f"{n}vecteurs_propres"), vecteur_propre)


def project_on_PCA_basis(datapath, n, v_min, v_max):
    """Project the raw data on the PCA basis.

    Args:
        datapath (str): Data path.
        n (int): Number of eigen vectors.
        v_min (float): Minimum value of the energy.
        v_max (float): Maximum value of the energy.
    """

    plot_spectra = True

    name_folder = f"{n}vecteurs_propres"
    basispath = os.path.join(datapath, name_folder)
    basisname = f"{n}vecteurs_propres.npy"

    # Get already corrected data
    name_data = "Corrected_spectra"
    name_wl = "Wavelength"
    stl_data = np.load(os.path.join(datapath, name_data + ".npy"))
    wavelength = np.load(os.path.join(datapath, name_wl + ".npy"))
    energy = 1240 / wavelength

    W1 = 1240 / v_min
    W2 = 1240 / v_max
    i1 = np.argmin(np.abs(wavelength - W1))
    i2 = np.argmin(np.abs(wavelength - W2))

    (m, l) = np.shape(stl_data)
    print(m, l)

    # Saving folder
    foldcontents = os.listdir(basispath)
    foldername = "Projection on PCA_basis"
    if foldername not in foldcontents:
        os.makedirs(os.path.join(basispath, foldername))
    savingpath = os.path.join(basispath, foldername)

    # Load and plot PCA basis
    PCABasis = np.load(os.path.join(basispath, basisname))
    lPCAB = np.shape(PCABasis)[1]

    couleur = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:brown",
        "tab:pink",
        "tab:gray",
        "tab:olive",
        "tab:cyan",
    ]

    for i in range(n):
        fig = plt.figure(figsize=[7, 5])
        ax = fig.add_subplot(111)

        plt.plot(energy[i1:i2], PCABasis[i][i1:i2], color=couleur[i % 1])

        plt.rcParams.update({"font.size": 8})
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=12)
        plt.xlabel("Energy (eV)", fontsize=14)
        plt.ylabel(f"PCA eigenvector #{i}", fontsize=14)

        plt.savefig(
            os.path.join(savingpath, f"vecteur_propre_#{i}.png"),
            dpi=300,
            bbox_inches="tight",
        )

    start = 0
    stop = m
    coeff_proj = []

    for j in range(start, stop, 1):
        print(j)

        spectrum = stl_data[j]

        if plot_spectra and j % 10 == 0:
            fig = plt.figure(figsize=[7, 5])
            ax = fig.add_subplot(111)

        fit = np.zeros(lPCAB)
        coeff = []
        s = ""

        for i in range(n):
            vector = PCABasis[i]
            a = np.dot(spectrum, vector)

            # if i == 2 or i==3 :
            #     a = 0

            if plot_spectra and j % 10 == 0:
                plt.plot(
                    energy[i1:i2],
                    a * vector[i1:i2],
                    label=f"f{i}",
                    color=couleur[i % 10],
                )

            fit += a * vector
            coeff.append(a * max(np.abs(vector)))
            s += f"c{i}={int(a * max(np.abs(vector)))}"  # round(a*max(vector), 2))

        coeff_proj.append(coeff)

        if plot_spectra and j % 10 == 0:
            s = s[:-2]
            plt.plot(energy[i1:i2], spectrum[i1:i2], "c", label="data")
            plt.plot(energy[i1:i2], fit[i1:i2], "k", label="Projection")

            plt.text(
                0.5,
                1.07,
                s,
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=12,
            )
            plt.rcParams.update({"font.size": 10})
            plt.legend(loc="upper right")
            plt.legend(bbox_to_anchor=(1, 1))
            ax.tick_params(axis="x", labelsize=16)
            ax.tick_params(axis="y", labelsize=16)
            plt.xlabel("Energy (eV)", fontsize=16)
            plt.ylabel("Intensity (arb. units)", fontsize=16)

            savename = os.path.join(savingpath, f"Spectrum{j}.png")
            plt.savefig(savename, dpi=300, bbox_inches="tight")

            plt.show()

    np.save(
        os.path.join(basispath, f"coeff_projection_from{start}to{stop}"), coeff_proj
    )
