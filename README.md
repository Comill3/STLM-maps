# STLM-maps
STLM acquisition data workflow

Authors: Mylène Sauty

Modifications: Camille Fornos

How the data should be saved on the computer:
.../STL_Data/XXXX-XX-XX Sample_name W UHV XX/STL_XX/...
Date, sample name, composition of the tip, UHV conditions, temperature (LN(liquide nitrogen) or RT(room temperature))

### 0_analyse_trigger.py ###
The spectrometer and the STM microscope communicate with triggers to acquire spectra. There is miscommunication between them leading to spectra not being acquired. They are called "missing spectra". Because the spectrometer didn't receive the trigger to stop the acquisition, the next spectrum will be "perverted" and be the result of averaging 2 pixels.

Identification of the missing spectra leads to the identification of the perverted spectra.

### 1_correct_and_save_spectra.py ###
Allows the identification and removing of cosmic rays in the spectra using wavelets. The program needs utils_EA_0.py, Clean_EA.py and the folder "filters" to operate.

### 2_plot_raw_spectra.py ###
Allow to plot the spectra without any noise treatement. Additionnaly, a plot of the summed spectra with its FWHM can be plotted.

### 3_remove_noise_PCA.py ###
Remove the noise by using the PCA method. Basically, we are creating a basis where we can project our raw data. Each vector of the basis is a composant of the sepctra and added to a specific constant, we can "fit" the data and have spectra with less noise.

### 4_plot_PCA_spectra.py ###
Allow to plot the spectra with noise treatement. Additionnaly, a plot of the summed spectra with its FWHM can be plotted.

### 5_create_gif.py ###
Create a gif to superpose on STM effective topographies. Each horizontal line is plotted and have their corresponding spectra.

### 6_fitting_data.py ###
Use of a combination of Gaussian function to fit the raw data. The position, width and amplitude of the Gaussain gives information on the QWs composition. Better to use on samples with one or two main emission peaks.

### 7_fits_maps.py ###
With the fitting, Gaussian coefficients are obtained (Amplitude, FWHW, energy position). They can be plotted spatially on a 2D map. This can be then associated with the effective STM topographies.

### 8_integrate_intensity.py ###
Integrate the PCA or raw data after the removing of cosmic rays. Then creates map with the corresponding histogram. Takes into account the missing spectra.


