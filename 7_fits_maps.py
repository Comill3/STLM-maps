"""Plot maps from the coefficient obtained thanks to the previous fitting."""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.cm as cm
from fit_functions import get_gaussian_coeff

datapath = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-03 #SG19 ALE W UHV LN/STL_37'

PCA = False
Log = False
n = 7
start = 0
stop = 1024
G = '8'

dataname, basispath = get_gaussian_coeff(PCA, datapath, n, G, start, stop)

# Get coeff Gaussian fit
coeff_Gaussian_fit = np.loadtxt(os.path.join(basispath, dataname))
(m,l) = np.shape(coeff_Gaussian_fit)
print(m, l)

#Saving folder
foldcontents = os.listdir(basispath)
foldername = 'Maps_Fit_Gaussian'
if foldername not in foldcontents:
    os.makedirs(os.path.join(basispath,foldername))
savingpath = os.path.join(basispath,foldername)
    
# Plot in grid
image_size = [500, 500] # [nm, nm]
grid_size = [32, 32] #32 #grid 
step_linescan = image_size[0]/grid_size[0] #in nm #for STL measurements
missing_spectra = []
perverted_spectra = [i+1 for i in missing_spectra] + [332, 333, 334, 335, 336, 337, 338, 339, 340]
    
# Add zeros for missing spectra
l2 = len(missing_spectra)
fit_param_modif = coeff_Gaussian_fit[:stop-start-l2, :]
zeroline = [0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
for i in missing_spectra :
    if i >= np.shape(fit_param_modif)[0]:
        fit_param_modif = np.append(fit_param_modif, np.array([zeroline]), axis=0)
    else :
        fit_param_modif = np.insert(fit_param_modif, i, zeroline, axis=0)
        
for i in perverted_spectra :
    fit_param_modif[i] = zeroline
        
print(np.shape(fit_param_modif))

### Quantum well emission Gaussian ###

A1 = np.reshape(np.copy(fit_param_modif[:, 1]), grid_size)
A1 = np.flip(np.abs(A1), axis = 1)
A1[A1<1] = 0

E1 = np.reshape(np.copy(fit_param_modif[:, 2]), grid_size)
E1 = np.flip(E1, axis = 1)

sigma1 = np.reshape(np.copy(fit_param_modif[:, 3]), grid_size)
sigma1 = np.flip(sigma1, axis = 1)

fig = plt.figure()
fig.set_figheight(3)
fig.set_figwidth(16)
spec = gridspec.GridSpec(ncols=3, nrows=1, width_ratios=[1, 1, 1], wspace=0.13)
cmap1 = cm.get_cmap(name='viridis', lut=None)
   
ax1 = fig.add_subplot(spec[0])

if Log:
    plt.imshow(np.log10(A1), cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
    clb = plt.colorbar(pad = 0.2)
    clb.set_label('log10(A1) (arb. units)', labelpad=20, y=0.5, rotation=270, fontsize=16)

else:
    plt.imshow(A1, cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
    clb = plt.colorbar(pad = 0.2)
    clb.set_label('A1 (arb. units)', labelpad=20, y=0.5, rotation=270, fontsize=16)

clb.ax.tick_params(labelsize=14) 
ax1.set_ylabel('Position (nm)', rotation=270, fontsize=12)
ax1.set_xlabel('Position (nm)', fontsize=12)
ax1.tick_params(axis='x', labelsize=10)
ax1.tick_params(axis='y', labelsize=10)
ax1.xaxis.tick_top()
ax1.yaxis.tick_right()
ax1.xaxis.set_label_coords(0.5, 1.25)
ax1.yaxis.set_label_coords(1.18, 0.5)

ax2 = fig.add_subplot(spec[1])
plt.imshow(E1, cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
clb = plt.colorbar(pad = 0.2) #, ticks=[2.30, 2.31, 2.32, 2.33, 2.34, 2.35])
clb.set_label('E1 (eV)', labelpad=20, y=0.5, rotation=270, fontsize=16)
clb.ax.tick_params(labelsize=14) 
ax2.set_ylabel('Position (nm)', rotation=270, fontsize=12)
ax2.set_xlabel('Position (nm)', fontsize=12)
ax2.tick_params(axis='x', labelsize=10)
ax2.tick_params(axis='y', labelsize=10)
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.xaxis.set_label_coords(0.5,1.25)
ax2.yaxis.set_label_coords(1.18,0.5)

ax3 = fig.add_subplot(spec[2])
plt.imshow(2.355*sigma1, cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
clb = plt.colorbar(pad = 0.2)
clb.set_label('$2.355 * \sigma$1 (eV)', labelpad=20, y=0.5, rotation=270, fontsize=16)
clb.ax.tick_params(labelsize=14) 
ax3.set_ylabel('Position (nm)', rotation=270, fontsize=12)
ax3.set_xlabel('Position (nm)', fontsize=12)
ax3.tick_params(axis='x', labelsize=10)
ax3.tick_params(axis='y', labelsize=10)
ax3.xaxis.tick_top()
ax3.yaxis.tick_right()
ax3.xaxis.set_label_coords(0.5,1.25)
ax3.yaxis.set_label_coords(1.18,0.5)   

savename = os.path.join(savingpath, G + 'G_Coeff_map_Gaussian1.png')
plt.savefig(savename, dpi = 300, bbox_inches='tight', format='png')
plt.show()

### Quantum well emission Gaussian - histogram ###
A1_h = np.copy(A1)
A1_h = A1_h.flatten()
E1_h = np.copy(E1)
E1_h = E1_h.flatten()
sigma1_h = np.copy(sigma1)
sigma1_h = sigma1_h.flatten()

fig = plt.figure()
fig.set_figheight(3)
fig.set_figwidth(15)
spec = gridspec.GridSpec(ncols=3, nrows=1, width_ratios=[1, 1, 1], wspace=0.25)

ax1 = fig.add_subplot(spec[0])
plt.hist(A1_h, bins=50)
ax1.set_ylabel('Histogram A1', fontsize=16)
ax1.set_xlabel('A1 (arb. units)', fontsize=16)
ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)

ax2 = fig.add_subplot(spec[1])
plt.hist(E1_h, bins=50)
ax2.set_ylabel('Histogram E1', fontsize=16)
ax2.set_xlabel('E1 (eV)', fontsize=16)
ax2.tick_params(axis='x', labelsize=14)
ax2.tick_params(axis='y', labelsize=14)

ax3 = fig.add_subplot(spec[2])
plt.hist(2.355 * sigma1_h, bins=50)
ax3.set_ylabel('Histogram 2.355 * $\sigma$1', fontsize=16)
ax3.set_xlabel('$2.355 * \sigma$1 (eV)', fontsize=16)
ax3.tick_params(axis='x', labelsize=14)
ax3.tick_params(axis='y', labelsize=14)

savename = os.path.join(savingpath, G + 'G_Hist_Gaussian1.png')
plt.savefig(savename, dpi = 300, bbox_inches='tight', format='png')
plt.show()
    
## Phonon replica emission Gaussian ###
        
# B1 = np.reshape(np.copy(fit_param_modif[:,4]), grid_size)
# B1 = np.flip(B1, axis = 1)
# B1[B1<1] = 0

# E1_Eph = np.reshape(np.copy(fit_param_modif[:,5]), grid_size)
# E1_Eph = np.flip(E1_Eph, axis = 1)

# sigma1 = np.reshape(np.copy(fit_param_modif[:,6]), grid_size)
# sigma1 = np.flip(sigma1, axis = 1)

# fig = plt.figure()
# fig.set_figheight(3)
# fig.set_figwidth(15)
# spec = gridspec.GridSpec(ncols=3, nrows=1, width_ratios=[1, 1, 1], wspace=0.13)
# cmap1 = cm.get_cmap(name='inferno', lut=None)

# ax1 = fig.add_subplot(spec[0])
# plt.imshow(np.log10(B1), cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
# clb = plt.colorbar(pad = 0.2)
# clb.set_label('log10(B1) (arb. units)', labelpad=20, y=0.5, rotation=270, fontsize=16)
# clb.ax.tick_params(labelsize=14)
# ax1.set_ylabel('Position (nm)', rotation=270, fontsize=12)
# ax1.set_xlabel('Position (nm)', fontsize=12)
# ax1.tick_params(axis='x', labelsize=10)
# ax1.tick_params(axis='y', labelsize=10)
# ax1.xaxis.tick_top()
# ax1.yaxis.tick_right()
# ax1.xaxis.set_label_coords(0.5,1.25)
# ax1.yaxis.set_label_coords(1.18,0.5)

# ax2 = fig.add_subplot(spec[1])
# plt.imshow(E1_Eph, cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
# clb = plt.colorbar(pad = 0.2)
# clb.set_label('E1 - Eph (eV)', labelpad=20, y=0.5, rotation=270, fontsize=16)
# clb.ax.tick_params(labelsize=14)
# ax2.set_ylabel('Position (nm)', rotation=270, fontsize=12)
# ax2.set_xlabel('Position (nm)', fontsize=12)
# ax2.tick_params(axis='x', labelsize=10)
# ax2.tick_params(axis='y', labelsize=10)
# ax2.xaxis.tick_top()
# ax2.yaxis.tick_right()
# ax2.xaxis.set_label_coords(0.5,1.25)
# ax2.yaxis.set_label_coords(1.18,0.5)

# ax3 = fig.add_subplot(spec[2])
# plt.imshow(2.355 * sigma1, cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
# clb = plt.colorbar(pad = 0.2)
# clb.set_label('$2.355 * \sigma$1 (eV)', labelpad=20, y=0.5, rotation=270, fontsize=16)
# clb.ax.tick_params(labelsize=14)
# ax3.set_ylabel('Position (nm)', rotation=270, fontsize=12)
# ax3.set_xlabel('Position (nm)', fontsize=12)
# ax3.tick_params(axis='x', labelsize=10)
# ax3.tick_params(axis='y', labelsize=10)
# ax3.xaxis.tick_top()
# ax3.yaxis.tick_right()
# ax3.xaxis.set_label_coords(0.5,1.25)
# ax3.yaxis.set_label_coords(1.18,0.5)

# savename = os.path.join(savingpath, G + 'G_Coeff_map_Gaussian1-phonon-replica.png')
# plt.savefig(savename, dpi = 300, bbox_inches='tight', format='png')
# plt.show()

### Phonon replica emission Gaussian - histogram ###

# B1_h = np.copy(B1)
# B1_h = B1_h.flatten()
# E1_Eph_h = np.copy(E1_Eph)
# E1_Eph_h = E1_Eph_h.flatten()
# sigma1_h = np.copy(sigma1)
# sigma1_h = sigma1_h.flatten()

# fig = plt.figure()
# fig.set_figheight(3)
# fig.set_figwidth(15)
# spec = gridspec.GridSpec(ncols=3, nrows=1, width_ratios=[1, 1, 1], wspace=0.25)

# ax1 = fig.add_subplot(spec[0])
# plt.hist(B1_h, bins=50)
# ax1.set_ylabel('Histogram B1', fontsize=16)
# ax1.set_xlabel('B1 (arb. units)', fontsize=16)
# ax1.tick_params(axis='x', labelsize=14)
# ax1.tick_params(axis='y', labelsize=14)

# ax2 = fig.add_subplot(spec[1])
# plt.hist(E1_Eph_h, bins=50)
# ax2.set_ylabel('Histogram E1 - Eph', fontsize=16)
# ax2.set_xlabel('E1 - Eph (eV)', fontsize=16)
# ax2.tick_params(axis='x', labelsize=14)
# ax2.tick_params(axis='y', labelsize=14)

# ax3 = fig.add_subplot(spec[2])
# plt.hist(2.355 * sigma1_h, bins=50)
# ax3.set_ylabel('Histogram 2.355 * $\sigma$1', fontsize=16)
# ax3.set_xlabel('$2.355 * \sigma$1 (eV)', fontsize=16)
# ax3.tick_params(axis='x', labelsize=14)
# ax3.tick_params(axis='y', labelsize=14)

# savename = os.path.join(savingpath, Gnumber + 'G_Hist_Gaussian1-phonon-replica.png')
# plt.savefig(savename, dpi = 300, bbox_inches='tight', format='png')
# plt.show()

# ### Ratio A1 to B1 ###

# B1 = np.reshape(np.copy(fit_param_modif[:, 4]), grid_size)
# B1 = np.flip(np.abs(B1), axis = 1)
# Ratio = B1/A1
# Ratio[A1==0] = 0

# fig = plt.figure()
# fig.set_figheight(3)
# fig.set_figwidth(5)
# spec = gridspec.GridSpec(ncols=1, nrows=1)
# cmap1 = cm.get_cmap(name='inferno', lut=None)
# cmap1.set_over('w', alpha=0.0)
# cmap1.set_under('w', alpha=0.0)

# ax1 = fig.add_subplot(spec[0])
# plt.imshow(np.abs(Ratio), cmap=cmap1, aspect=1, extent=[image_size[1], 0, image_size[0], 0], interpolation='none')
# clb = plt.colorbar(pad = 0.2)
# clb.set_label('B1/A1 (arb. units)', labelpad=20, y=0.5, rotation=270, fontsize=16)
# clb.ax.tick_params(labelsize=14)
# ax1.set_ylabel('Position (nm)', rotation=270, fontsize=12)
# ax1.set_xlabel('Position (nm)', fontsize=12)
# ax1.tick_params(axis='x', labelsize=10)
# ax1.tick_params(axis='y', labelsize=10)
# ax1.xaxis.tick_top()
# ax1.yaxis.tick_right()
# ax1.xaxis.set_label_coords(0.5,1.25)
# ax1.yaxis.set_label_coords(1.18,0.5)

# savename = os.path.join(savingpath, Gnumber + 'G_Coeff_map_B1divA1.png')
# plt.savefig(savename, dpi = 300, bbox_inches='tight', format='png')
# plt.show()

# ### Ratio A2 to A1 - histogram ###
# B1sA1_h = np.copy(B1/A1)
# B1sA1_h[B1sA1_h > 1] = np.nan
# B1sA1_h = B1sA1_h.flatten()

# fig = plt.figure()
# fig.set_figheight(3)
# fig.set_figwidth(5)
# spec = gridspec.GridSpec(ncols=1, nrows=1)

# ax1 = fig.add_subplot(spec[0])
# plt.hist(B1sA1_h, bins=50)
# ax1.set_ylabel('Histogram B1/A1', fontsize=16)
# ax1.set_xlabel('B1/A1 (arb. units)', fontsize=16)
# ax1.tick_params(axis='x', labelsize=14)
# ax1.tick_params(axis='y', labelsize=14)

# savename = os.path.join(savingpath, Gnumber + 'G_Hist_B1divA1.png')
# plt.savefig(savename, dpi = 300, bbox_inches='tight', format='png')
# plt.show()