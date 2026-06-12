"""Removes noise from the raw data using the Principal Componants Ananlysis method."""

from PCA_functions import perform_PCA

datapath = "Demo_Data/STL_8"
n = 7  # number of eigen vectors
v_min = 1.5
v_max = 3.5

### Uncomment first until the right number of eigne vector is found ###

# Find the number of eigen vectors
perform_PCA(datapath, n)

### Uncomment second ###

# # Projection of the data in the basis
# project_on_PCA_basis(datapath, n, v_min, v_max)
