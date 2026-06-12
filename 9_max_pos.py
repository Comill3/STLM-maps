"""Plot maps from the coefficient obtained thanks to the previous fitting."""

import os

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

from fit_functions import get_data_for_max

datapath = "Demo_Data/STL_8"

PCA = True
n = 7
start = 0
stop = 1024

# mpl.rcParams['font.family'] = 'Arial'

if PCA:
    save_path = os.path.join(datapath, f"{n}vecteurs_propres")
else:
    save_path = datapath

# Saving folder
foldcontents = os.listdir(save_path)
foldername = "Maps_max_pos"
if foldername not in foldcontents:
    os.makedirs(os.path.join(save_path, foldername))
savingpath = os.path.join(save_path, foldername)

# Get the data
stl_data, stl_data_init, wavelength, energy = get_data_for_max(PCA, datapath, n)

# Plot in grid
image_size = [2000, 2000]  # [nm, nm]
grid_size = [32, 32]  # 32 #grid
step_linescan = image_size[0] / grid_size[0]  # in nm #for STL measurements
missing_spectra = [55, 260, 412, 429, 612]
perverted_spectra = [i + 1 for i in missing_spectra]

max_pos = []

W1 = 1240 / 1.5
W2 = 1240 / 3.5
i1 = np.argmin(np.abs(wavelength - W1))
i2 = np.argmin(np.abs(wavelength - W2))
#    i1 = 0
#    i2 = -1

liste_spectra = np.arange(0, 1024, 1)

for j in liste_spectra:
    print(j)

    spectrum = stl_data[j]
    spectrum = savgol_filter(spectrum, 20, 1)

    spectrum_init = stl_data_init[j]

    # Find maximum position

    A0 = np.max(spectrum[i1:i2])
    E0 = energy[i1 + np.argmax(spectrum[i1:i2])]

    max_pos.append([A0, E0])

max_pos = np.array(max_pos)

(m, l) = np.shape(max_pos)
print(m, l)

# Add zeros for missing spectra
l2 = len(missing_spectra)
max_pos_modif = max_pos[: stop - start - l2, :]
zeroline = [np.nan, np.nan]
for i in missing_spectra:
    if i >= np.shape(max_pos_modif)[0]:
        max_pos_modif = np.append(max_pos_modif, np.array([zeroline]), axis=0)
    else:
        max_pos_modif = np.insert(max_pos_modif, i, zeroline, axis=0)

for i in perverted_spectra:
    max_pos_modif[i] = zeroline

print(np.shape(max_pos_modif))

### Quantum well emission Gaussian ###

A0 = np.reshape(np.copy(max_pos_modif[:, 0]), grid_size)
A0 = np.flip(np.abs(A0), axis=1)

E0 = np.reshape(np.copy(max_pos_modif[:, 1]), grid_size)
E0 = np.flip(E0, axis=1)

fig, ax = plt.subplots()

cmap1 = cm.get_cmap(name="viridis", lut=None)

plt.imshow(
    E0,
    cmap=cmap1,
    aspect=1,
    extent=[image_size[1], 0, image_size[0], 0],
    interpolation="none",
)
clb = plt.colorbar(pad=0.2)
clb.set_label("E0 (eV)", labelpad=20, y=0.5, rotation=270, fontsize=16)

clb.ax.tick_params(labelsize=14)
ax.set_ylabel("Position (nm)", rotation=270, fontsize=14)
ax.set_xlabel("Position (nm)", fontsize=14)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
ax.xaxis.tick_top()
ax.yaxis.tick_right()
ax.xaxis.set_label_coords(0.5, 1.25)
ax.yaxis.set_label_coords(1.18, 0.5)

savename = os.path.join(savingpath, f"Max_energy_pos_{energy[i1]}_{energy[i2]}eV.png")
plt.savefig(savename, dpi=300, bbox_inches="tight", format="png")
plt.show()

fig, ax = plt.subplots()

cmap1 = cm.get_cmap(name="viridis", lut=None)

plt.imshow(
    np.log10(A0),
    cmap=cmap1,
    aspect=1,
    extent=[image_size[1], 0, image_size[0], 0],
    interpolation="none",
)
clb = plt.colorbar(pad=0.2)
clb.set_label("log10(A0)", labelpad=20, y=0.5, rotation=270, fontsize=16)

clb.ax.tick_params(labelsize=14)
ax.set_ylabel("Position (nm)", rotation=270, fontsize=14)
ax.set_xlabel("Position (nm)", fontsize=14)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
ax.xaxis.tick_top()
ax.yaxis.tick_right()
ax.xaxis.set_label_coords(0.5, 1.25)
ax.yaxis.set_label_coords(1.18, 0.5)

savename = os.path.join(savingpath, f"Max_amp_{energy[i1]}_{energy[i2]}eV.png")
plt.savefig(savename, dpi=300, bbox_inches="tight", format="png")
plt.show()

### Quantum well emission Gaussian - histogram ###
E0_h = np.copy(E0)
E0_h = E0_h.flatten()

A0_h = np.copy(A0)
A0_h = A0_h.flatten()

fig, ax = plt.subplots()

plt.hist(E0_h, bins=50)
ax.set_ylabel("Histogram E0", fontsize=16)
ax.set_xlabel("E0 (eV)", fontsize=16)
ax.tick_params(axis="x", labelsize=14)
ax.tick_params(axis="y", labelsize=14)

fig, ax = plt.subplots()

plt.hist(A0_h, bins=50)
ax.set_ylabel("Histogram A0", fontsize=16)
ax.set_xlabel("A0 (arb. units)", fontsize=16)
ax.tick_params(axis="x", labelsize=14)
ax.tick_params(axis="y", labelsize=14)

# savename = os.path.join(savingpath, G + 'G_Hist_Gaussian1.png')
# plt.savefig(savename, dpi = 300, bbox_inches='tight', format='png')
# plt.show()
