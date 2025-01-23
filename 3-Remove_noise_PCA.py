"""Removes noise from the raw data using the Principal Componants Ananlysis method."""

from PCA_functions import perform_PCA, project_on_PCA_basis

datapath = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-03 #SG19 ALE W UHV LN/STL_37'
n = 7
min = 1.5
max = 3.5

# Find the number of eigen vectors
perform_PCA(datapath, n)

# # Projection of the data in the basis
# project_on_PCA_basis(datapath, n, min, max)
