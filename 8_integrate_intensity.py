"""Integrate pectra intensity from raw or PCA data."""

import os

import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from plot_functions import integrate_data

# mpl.rcParams['font.family'] = 'Arial'

DATA_PATH = "Demo_Data"
DEDICATED_FOLDER = "STL_8"
DATA_PATH2 = os.path.join(DATA_PATH, DEDICATED_FOLDER)

PCA = True
Smooth = False
n = 7  # number of eigen vectors

# Intg limits
E1 = 3.01  # eV
E2 = 3.51  # eV

# Acquisition parameters #
image_size = [1000, 1000]
grid_size = [32, 32]
missing_spectra = [0]
perverted_spectra = [i + 1 for i in missing_spectra]

# ## Integrate data and save it in a text file ###
integrate_data(DATA_PATH2, PCA, n, Smooth, E1, E2, grid_size[0], grid_size[1])

### Plot and save the corresponding map ###
Ech_log = True
start = 0
stop = 1024

file_name = f"Intg_matrix_from_{E1}_to_{E2}eV.txt"

if PCA:
    name_folder = f"{n}vecteurs_propres"
    data_path = os.path.join(DATA_PATH2, name_folder)
else:
    data_path = DATA_PATH2

data_path2 = os.path.join(data_path, "Correlations")

data_intg = np.loadtxt(os.path.join(data_path2, file_name))
data_intg = np.flip(data_intg, axis=1)
data_intg = data_intg.flatten().tolist()

l2 = len(missing_spectra)
fit_param_modif = data_intg[: stop - start - l2]
zeroline = [np.nan]
for i in missing_spectra:
    if i >= np.shape(fit_param_modif)[0]:
        fit_param_modif = np.append(fit_param_modif, np.array([zeroline]), axis=0)
    else:
        fit_param_modif = np.insert(fit_param_modif, i, zeroline, axis=0)

for i in perverted_spectra:
    fit_param_modif[i] = np.nan

I = np.reshape(np.copy(fit_param_modif), grid_size)
I = np.flip(np.abs(I), axis=1)

if Ech_log:
    I[I < 1] = 1

fig = plt.figure()
spec = gridspec.GridSpec(ncols=1, nrows=1)
cmap1 = cm.get_cmap(name="viridis", lut=None)

ax1 = fig.add_subplot(spec[0])
if Ech_log:
    plt.imshow(
        np.log10(I),
        cmap=cmap1,
        extent=[image_size[1], 0, image_size[0], 0],
        interpolation="none",
    )
    clb = plt.colorbar(pad=0.2)
    clb.set_label(
        "log10(Integrated intensity)", labelpad=20, y=0.5, rotation=270, fontsize=16
    )

else:
    plt.imshow(
        I, cmap=cmap1, extent=[image_size[1], 0, image_size[0], 0], interpolation="none"
    )
    clb = plt.colorbar(pad=0.2)
    clb.set_label(
        "Integrated intensity (arb. units)",
        labelpad=20,
        y=0.5,
        rotation=270,
        fontsize=16,
    )

# clb.ax.tick_params(labelsize=14)
ax1.set_ylabel("Position (nm)", rotation=270, fontsize=16)
ax1.set_xlabel("Position (nm)", fontsize=16)
ax1.tick_params(axis="x", labelsize=12)
ax1.tick_params(axis="y", labelsize=12)
ax1.xaxis.tick_top()
ax1.yaxis.tick_right()
ax1.xaxis.set_label_coords(0.5, 1.25)
ax1.yaxis.set_label_coords(1.18, 0.5)

if Ech_log:
    save_name = file_name[:-4] + "_log"
else:
    save_name = file_name[:-4]

savename = os.path.join(data_path2, save_name + ".png")
plt.savefig(savename, dpi=300, bbox_inches="tight", format="png")
plt.show()

### Histogram ###
I_h = np.copy(I)
I_h = I_h.flatten()

fig = plt.figure()

plt.hist(I_h, bins=100)
plt.ylabel("Histogram", fontsize=16)

plt.xlabel("I (arb. units)", fontsize=16)

plt.tick_params(axis="x", labelsize=14)
plt.tick_params(axis="y", labelsize=14)
