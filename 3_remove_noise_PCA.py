"""Removes noise from the raw data using the Principal Componants Ananlysis method."""

from PCA_functions import perform_PCA, project_on_PCA_basis

datapath = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/' \
            '2024-07-18 #SG19 ALE W UHV LN/STL_14'
n = 10
v_min = 1.5
v_max = 3.5

# # Find the number of eigen vectors
# perform_PCA(datapath, n)

# Projection of the data in the basis
project_on_PCA_basis(datapath, n, v_min, v_max)
