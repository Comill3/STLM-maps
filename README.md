# STLM-maps
STLM acquisition data workflow

### 0_analyse_trigger.py ###
The spectrometer and the STM microscope communicate with triggers to acquire spectra. There is miscommunication between them leading to spectra not being acquired. They are called "missing spectra". Because the spectrometer didn't receive the trigger to stop the acquisition, the next spectrum will be "perverted" and be the result of averaging 2 pixels.

Identification of the missing spectra leads to the identification of the perverted spectra.

### 1_correct_and_save_spectra.py ###
Allows the identification and removing of cosmic rays in the spectra using wavelets. The program neads utils_EA_0.py and the folder "filters" to operate.

### 2_plot_raw_spectra.py ###
