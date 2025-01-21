#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 21:27:51 2021

@author: berlioz
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def LoadWavelet1D(j, N):
    """
    Generate the fourier transform of a bump steerable wavelet at a 2^j
    characteristic scale, defined on an array of size N
    """
    
    # Chargement de la banque d'ondeltte pré-calculée
    J, L, dn, df = 2, 4, 0,  3 #2, 4, 0, 3 #3, 8, 0,  3
    filename = 'bump_steerable_wavelet' + '_N' + str(N) + '_J' + str(J) + '_L' + str(L) + '_dn' + str(dn) + '_df' + str(df)  + '.npy'
    matfilters = np.load(os.path.join('filters',filename), allow_pickle=True).item()
    fftpsi = matfilters['filt_fftpsi'].astype(np.complex_)
    
    # Selection de l'ondelette
    l1 = 0
    psi_hat = fftpsi[l1,j,0,:,:]
    
    # Passage de 2d à 1d
    wave1 = np.fft.ifft2(psi_hat)
    wave1 = wave1[0,:]
    wave1 = np.fft.fft(wave1)
    
#    wavelet_to_plot = np.fft.ifft(wave1)
#    wavelet_to_plot = [wavelet_to_plot[i-20] for i in range(len(wavelet_to_plot))]
#    fig = plt.figure(1)
#    ax = fig.add_subplot(111)
#    plt.plot(wavelet_to_plot, 'r')
#    ax.tick_params(axis='x', labelsize=16)
#    ax.tick_params(axis='y', labelsize=16)
#    plt.xlim(0,40)
#    plt.show()
    
    return wave1

def Symmetrize1D(data_ini, N):
    """
    Symmetrize a 1D array of size < N into an array of size N, by the right
    """
    
    # On calcul le nombre de fois qu'il faudra flipper le vecteur initial
    shape = data_ini.shape[0]
    nb_flip = int(np.ceil(N/shape))
    
    # On remplit le vecteur de taille nb_flip*shape>N en faisant des flips successifs
    data = np.zeros((nb_flip*shape))
    odd = 1
    for i in range(nb_flip):
        if odd == 1:
            data[i*shape:(i+1)*shape] = data_ini
            odd = 0
        elif odd == 0:
            data[i*shape:(i+1)*shape] = np.flip(data_ini, axis=0)
            odd = 1
    
    # On recoupe en taille N
    data = data[:N]
    
    return data

def CleanLine(data_line, wavelet, eps, seuil_bas):
    """
    Clean a line and return the position of all cosmic rays
    Need a 1D wavelet, the line will be extended to the size
    of the wavelet. The size of the wavelet has to be larger 
    than the size of the line.
    """
#    fig = plt.figure(2)
#    ax = fig.add_subplot(111)
#    plt.plot(data_line)
#    ax.tick_params(axis='x', labelsize=16)
#    ax.tick_params(axis='y', labelsize=16)
#    plt.show()
    
    N = wavelet.shape[0]
    size_ini = data_line.shape[0]
    
    error = False
        
    PosRC = []
    crit = 1
    step = 10
    while crit == 1 and step > 0 :
        shape = data_line.shape[0] # on sauve la nouvelle taille
        data_line = Symmetrize1D(data_line, N) # on symetrise
        Conv1 = np.abs(np.real(np.fft.ifft(np.fft.fft(np.log(data_line))*wavelet))) # on convolue + abs
        Conv1 = Conv1[2:shape] # on recoupe en excluant les premiers points
        
#        fig = plt.figure(3)
#        ax = fig.add_subplot(111)
#        plt.plot(Conv1, 'k')
#        ax.tick_params(axis='x', labelsize=16)
#        ax.tick_params(axis='y', labelsize=16)
#        plt.show()
        
        maxval, argmax = np.amax(Conv1), np.argmax(Conv1) # on trouve le max et l'argmax
        argmax = argmax + 2  # on rajoute le 2 qu'on a coupé plus tôt
        data_line = data_line[:shape] # on recoupe les données à taille initiale
        if maxval > eps:
            data_loc = data_line[argmax-1: argmax+2] # on fait une étude locale
            # on commence par tester si le point le plus bas est en dessous du seuil bas
            minval, delta_arg = np.amin(data_loc), np.argmin(data_loc) # on trouve le min local et son arg
            # si minval est en dessous du seuil bas, on garde delta_arg
            # sinon on prend celui du max
            if minval > seuil_bas:
                delta_arg = np.argmax(data_line[argmax-1: argmax+2]) # on cherche localement la valeur max
            data_line = np.delete(data_line, argmax + delta_arg - 1) # on enleve le point en question
            PosRC.append(argmax + delta_arg - 1) # on enregistre sa position
        else:
            crit = 0
        step = step - 1
        if step == 0:
            error = True
            
    # On transforme PosRC pour correspondre au vecteur non-retranché 
    PosRC = PosProgToPosIni(PosRC, size_ini) 
    
    return PosRC, data_line, error

def PosProgToPosIni(PosRC, size_ini):
    """ 
    On obtient le vecteur de position des différents rayons cosmiques sur la line 
    initiale à partir de positions extraites progressivement. On reconstruit pour
    ça le vecteur, à défaut d'une stratégie plus maline...
    """
    
    fakedata = np.zeros(size_ini-len(PosRC))
    for i in range(len(PosRC)):
        fakedata = np.insert(fakedata, PosRC[-i-1], 1)
    PosRC = np.where(fakedata[:] == 1)
        
    return list(PosRC[0])


