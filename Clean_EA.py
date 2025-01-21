#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 21:10:27 2021

@author: berlioz
"""

import matplotlib.pyplot as plt
import numpy as np
import os
from utils_EA_0 import LoadWavelet1D, CleanLine

""" On charge les données """

# On charge les données
data_ini = np.load(os.path.join("Data","data_MLPB.npy"))

# On les plot
fig, ax = plt.subplots(1,1)
img = ax.imshow(np.log(data_ini), cmap = 'jet')
fig.colorbar(img)

""" Parametres et chargement des ondelettes """

# parametres
eps = 0.03
seuil_bas = 620
j = 0
N = 256

# chargement ondelette
wave1 = LoadWavelet1D(j, N)

""" Nettoyage d'une ligne """

# Selection de la ligne
#data_line_ini = data_ini[483,:]
data_line_ini = data_ini[342,:]

# on plot pour dev 
data_line = data_line_ini
plt.figure(1)
plt.plot(data_line)

PosRC, data_clean, error = CleanLine(data_line, wave1, eps, seuil_bas)

plt.figure(2)
plt.plot(data_clean)

""" On effectue le nettoyage sur l'ensemble de l'image """

data_clean_full = np.zeros((data_ini.shape[0],data_ini.shape[1]))
for i in range(data_ini.shape[0]):
    data_line = data_ini[i,:]
    _, data_clean, error = CleanLine(data_line, wave1, eps, seuil_bas)
    data_clean_full[i, :data_clean.shape[0]] = data_clean
    if error:
        print("something bad happened at line " + str(i))
    
# plot pour visualisation
plt.figure(1)
for i in range(data_ini.shape[0]):
    plt.plot(data_ini[i])

plt.figure(2)
for i in range(data_ini.shape[0]):
    plt.plot(data_clean_full[i])
