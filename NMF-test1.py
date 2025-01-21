# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:16:49 2022

@author: Mylène Sauty
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF




def Perform_NMF(datapath, Plot_spectra):    
    name_data = 'Selected_spectra_for_PCA'
    name_wl = 'Wavelength'

        
    Selected_spectra = np.load(os.path.join(datapath, name_data + '.npy'))
    print(np.shape(Selected_spectra))
    Wavelength = np.load(os.path.join(datapath, name_wl + '.npy'))
    Energy = 1240/Wavelength
    
    STL_data = Selected_spectra
    
    
    #Get data projected on PCA basis
    n = 6
    name_folder = '{}vecteurs_propres'.format(n)
    basispath = os.path.join(datapath,name_folder)
    basisname = "{}vecteurs_propres.npy".format(n)
    dataname = 'coeff_projection_from0to{}.npy'.format(1024)
    coeff = np.load(os.path.join(basispath, dataname))
    vectors = np.load(os.path.join(basispath, basisname))
    for i in range(np.shape(vectors)[0]):
        vectors[i] = vectors[i]/max(np.abs(vectors[i]))
    STL_data = np.matmul(coeff, vectors)
    print(np.shape(STL_data))
    
    (m, l) = np.shape(STL_data)    


    #Saving folder
    foldcontents = os.listdir(basispath)
    foldername = 'Projection on NMF_vectors'
    if foldername not in foldcontents:
        os.makedirs(os.path.join(basispath,foldername))
    savingpath = os.path.join(basispath,foldername)
    
    


    STL_data[STL_data<0] = 0

    nn = 8
    
    model = NMF(n_components=nn, init='random', random_state=0)
    W = model.fit_transform(STL_data)
    H = model.components_
    
    print(np.shape(W))
    print(np.shape(H))
    
    fig = plt.figure(figsize = [7, 5])
    ax = fig.add_subplot(111)
    plt.rcParams.update({'font.size': 8})
    for j in range(nn) :
        plt.plot(Energy, H[j], label=f'f{j}')
    plt.xlim(1.5, 3.5)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    plt.xlabel('Energy (eV)', fontsize=14)
    plt.ylabel('NMF vectors'.format(i), fontsize=14)
    plt.legend(loc='upper right')
    plt.legend(bbox_to_anchor=(1, 1))
    plt.savefig(os.path.join(savingpath,'NMF_vectors.png'), dpi = 300, bbox_inches='tight')
    plt.show()
    
    
    plt.plot(STL_data.T)
    #plt.xlim(700, 900)
    plt.savefig(os.path.join(savingpath,'test.png'), dpi = 300, bbox_inches='tight')
    plt.show()
    
    coeff_proj = []
    start = 0
    stop = m
    
    for i in range(start, stop) :
        coeff = []
        fit = 0
        s = ''
        if Plot_spectra :
            fig = plt.figure(figsize = [5, 4])
            plt.rcParams.update({'font.size': 12})
            plt.plot(Energy, STL_data[i], 'r')
        for j in range(nn) :
            vect = H[j]
            a = W[i,j]
            if Plot_spectra :
                plt.plot(Energy, a*vect, label='f{}'.format(j))
                fit += a*vect
                s += 'c{}={}; '.format(j, round(max(a*vect)))
            coeff.append(max(a*vect))
            
        if Plot_spectra :
            plt.plot(Energy, fit, 'k')
            plt.xlim(1.5, 3.5)
            plt.text(0.5, 1.2, s, ha='center', va='center', transform=ax.transAxes, fontsize=12)
            plt.legend(loc='upper right')
            plt.legend(bbox_to_anchor=(1, 1))
            ax.tick_params(axis='x', labelsize=14)
            ax.tick_params(axis='y', labelsize=14)
            plt.xlabel('Energy (eV)', fontsize=14)
            plt.ylabel('Intensity (arb. units)', fontsize=14)
            savename = os.path.join(savingpath,'Spectrum{}.png'.format(i))
            plt.savefig(savename, dpi = 100, bbox_inches='tight')
            plt.show()
        
        coeff.append(max(STL_data[i]))
        
        coeff_proj.append(coeff)
        
    np.save(os.path.join(basispath, 'coeff_projection_NMF_from{}to{}'.format(start, stop)), coeff_proj)

    
    
    

datapath = "C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-24 #SG19 ALE W UHV LN/STL_8"
Plot_spectra = True

Perform_NMF(datapath, Plot_spectra)


