"""Plot function."""

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.signal import savgol_filter
from scipy import integrate as intg

# mpl.rcParams['font.family'] = 'Arial'

def param_plot_2d(i, ax, xlabel, ylabel, s, date, spectrum_number, savename):
    """Plots a 2D plot with specified labels and saves the figure.

    Args:
        i (int): Figure number.
        ax (matplotlib.axes.Axes): Axes object to plot on.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        s (str): Text to display at the top center of the plot.
        date (str): Date to display at the top left of the plot.
        spectrum_number (int): Spectrum number to display next to the date.
        savename (str): Path and filename where the plot will be saved.
    """

    plt.figure(i)

    plt.text(0.75, 1.07, s, ha='center', va='center',
             transform=ax.transAxes, fontsize=12)
    plt.text(0.10, 1.07, str(date) + ' #' + str(spectrum_number),
             ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.savefig(savename, dpi=300, bbox_inches='tight')

def param_plot_3d(i, ax, xlabel, ylabel, zlabel, s, date, spectrum_number, savename):
    """Plots a 3D plot with specified labels and saves the figure.

    Args:
        i (int): Figure number.
        ax (matplotlib.axes.Axes): Axes object to plot on.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        zlabel (str): Label for the z-axis.
        s (str): Text to display at the top center of the plot.
        date (str): Date to display at the top left of the plot.
        spectrum_number (int): Spectrum number to display next to the date.
        savename (str): Path and filename where the plot will be saved.
    """
    plt.figure(i)
    ax.set_zlabel(zlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=12)

    # Change the backgroung color
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    # Additional text on the spectra
    # Place text on a fixed position on the ax object
    ax.text2D(0.4, 0.95, s, transform=ax.transAxes, fontsize=10)
    ax.text2D(0.15, 0.95, str(date) + ' #' + str(spectrum_number),
              transform=ax.transAxes, fontsize=10)

    plt.rcParams.update({'font.size': 8}) # Set the font size of the text
    # Set the elevation and azimuth of the axes in degrees (not radians)
    ax.view_init(35, -80)
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='z', labelsize=10)
    plt.savefig(savename, dpi=300, bbox_inches='tight')

def plot_stl_in_grid(stl_to_plot, datapath2, param_dict, step_linescan, grid_size, missing_spectra, perverted_spectra, Full, v_min, v_max, Smooth, SVG):
    """Plots STL spectra in a grid of specified size, with options for smoothing, 
    handling missing and perverted spectra, and plotting in both 2D and 3D.
    
    Args:
        stl_to_plot (str): The name of the STL file to plot.
        datapath2 (str): The path to the directory containing the data files.
        param_dict (dict): Dictionary containing parameters for each spectrum.
        step_linescan (float): The step size for the linescan.
        grid_size (int): The size of the grid for plotting.
        missing_spectra (list): List of indices of missing spectra.
        perverted_spectra (list): List of indices of perverted spectra.
        Full (bool): Whether to plot the full spectrum or a subset.
        v_min (float): Minimum value for wavelength or energy range.
        v_max (float): Maximum value for wavelength or energy range.
        Smooth (bool): Whether to apply smoothing to the spectra.
        SVG (bool): If True, save figs in .svg format, else save in .png format
    """

    date = stl_to_plot[:6]
    spectrum_number = int(stl_to_plot[7:-8])
    print(spectrum_number)

    with open(os.path.join(datapath2, stl_to_plot), 'r', encoding='utf-8') as entete:
        int_time = int(entete.readlines()[5].split()[17][:-1])
    
    # entete = open(os.path.join(datapath2, stl_to_plot), 'r')
    # int_time = int(entete.readlines()[5].split()[17][:-1])

    # Load data with cosmic rays already removed
    name_data = 'Corrected_spectra' #'Corrected_spectra' #'Interference_Corrected_spectra'
    name_wl = 'Wavelength'

    stl_data = np.load(os.path.join(datapath2, name_data + '.npy'))
    wavelength = np.load(os.path.join(datapath2, name_wl + '.npy'))
    (m, l) = np.shape(stl_data)

    # Add lines of zeros for missing spectra
    l2 = len(missing_spectra)
    spectra_modif = stl_data[0:m-l2]
    zeroline = np.zeros(l)
    for i in missing_spectra:
        if i > m-l2-1:
            spectra_modif = np.append(spectra_modif, [zeroline], axis=0)
        else:
            spectra_modif = np.insert(spectra_modif, i, zeroline, axis=0)

    # Reshape plots
    if not Full:
        W1 = 1240/v_min
        W2 = 1240/v_max
        i1 = np.argmin(np.abs(wavelength-W1))
        i2 = np.argmin(np.abs(wavelength-W2))

    else:
        i1 = 0
        i2 = -1

    # Plot STL spectra in grid of size: grid_size[0]*grid_size[1]

    for g in range(grid_size):

        num1 = g*grid_size
        num2 = (g+1)*grid_size

        energy = 1240/wavelength
        spectra = spectra_modif[num1:num2]

    # Plot the data: create the 4 figures for plot

        fig2 = plt.figure(2, figsize=[9, 6])
        ax2 = fig2.add_subplot(111, projection='3d')

        fig3 = plt.figure(3, figsize=[9, 6])
        ax3 = fig3.add_subplot(111, projection='3d')

        fig4 = plt.figure(4, figsize=[7, 5])
        ax4 = fig4.add_subplot(111)

        fig5 = plt.figure(5, figsize=[7, 5])
        ax5 = fig5.add_subplot(111)

        for k in range(num2-1, num1-1, -1):

            # remove perverted spectra and lines of zeros added at missing spectra places
            if k in missing_spectra or k in perverted_spectra:
                continue

            j = k - num1
            couleur = j/(num2-num1)
            # cmap = plt.cm.jet
            cmap = plt.cm.rainbow

            # Take one spectrum and remove constant background is needed
            intj = spectra[j][:]

            if Smooth:
                intj = savgol_filter(intj, 5, 1)

        # Create third direction for 3D graph
            Y = np.empty(len(intj))
            Y.fill((j+0.5)*step_linescan)

            # Corrected data vs wl 3D
            plt.figure(2)
            plt.plot(wavelength[i1:i2], Y[i1:i2], intj[i1:i2],
                     '-', color=cmap(couleur), linewidth=0.7)

            # Corrected data vs En 3D
            plt.figure(3)
            plt.plot(energy[i1:i2], Y[i1:i2], intj[i1:i2],
                     '-', color=cmap(couleur), linewidth=0.7)

            # Corrected data vs wl 2D
            plt.figure(4)
            plt.plot(wavelength[i1:i2], intj[i1:i2], '-',
                     color=cmap(couleur), linewidth=1, label='#{}'.format(Y[0]))

            # Corrected data vs En 2D
            plt.figure(5)
            plt.plot(energy[i1:i2], intj[i1:i2], '-',
                     color=cmap(couleur), linewidth=1, label='#{}'.format(Y[0]))

        # Plotting parameters for graphs and saving
        if spectrum_number in param_dict:
            s = f'line{g} ; {param_dict[spectrum_number][0]}V ; {param_dict[spectrum_number][1]}nA ; VLED = {param_dict[spectrum_number][2]}V ; {int_time}ms'
        else:
            s = f'{int_time}ms'

        # Folder name for saving plots
        foldcontents = os.listdir(datapath2)
        foldername = 'Plot_spectra'
        if Smooth:
            foldername += '_Smooth'

        if foldername not in foldcontents:
            os.makedirs(os.path.join(datapath2, foldername))
        savepath = os.path.join(datapath2, foldername)

        i = 2
        ax = ax2
        xlabel = 'Wavelength (nm)'
        ylabel = 'Position (nm)'
        zlabel = 'Intensity (arb. units)'
        savename = os.path.join(savepath, stl_to_plot.rstrip(
            '.txt') + f'vsWl_waterfall_from{num1}to{num2}')

        param_plot_3d(i, ax, xlabel, ylabel, zlabel, s,
                      date, spectrum_number, savename)

        i = 3
        ax = ax3
        xlabel = 'Energy (eV)'
        ylabel = 'Position (nm)'
        zlabel = 'Intensity (arb. units)'
        savename = os.path.join(savepath, stl_to_plot.rstrip(
            '.txt') + f'vsEn_waterfall_from{num1}to{num2}')

        param_plot_3d(i, ax, xlabel, ylabel, zlabel, s,
                      date, spectrum_number, savename)

        i = 4
        ax = ax4
        xlabel = 'Wavelength (nm)'
        ylabel = 'Intensity (arb. units)'
        savename = os.path.join(savepath, stl_to_plot.rstrip(
            '.txt') + f'vsWl_2D_from{num1}to{num2}')

        param_plot_2d(i, ax, xlabel, ylabel, s, date,
                      spectrum_number, savename)

        i = 5
        ax = ax5
        xlabel = 'Energy (eV)'
        ylabel = 'Intensity (arb. units)'

        if SVG:
            savename = os.path.join(savepath, stl_to_plot.rstrip(
                '.txt') + f'vsEn_2D_from{num1}to{num2}.svg')
        else:
            savename = os.path.join(savepath, stl_to_plot.rstrip(
                '.txt') + f'vsEn_2D_from{num1}to{num2}')

        param_plot_2d(i, ax, xlabel, ylabel, s, date,
                      spectrum_number, savename)

        plt.close('all')

def plot_one_spectrum(X, Y, couleur, xlabel, ylabel, s, s0, savename, ax):
    """Plots a single spectrum and saves the figure.
    
    Args:
        X (array-like): The x-axis data.
        Y (array-like): The y-axis data.
        couleur (str): The color of the plot line.
        xlabel (str): The label for the x-axis.
        ylabel (str): The label for the y-axis.
        s (str): Text to display at a specific location on the plot.
        s0 (str): Text to display at another specific location on the plot, in red color.
        savename (str): The filename to save the plot.
        ax (matplotlib.axes.Axes): The axes on which to plot.
    """

    plt.plot(X, Y, couleur)
    plt.xlabel(xlabel, fontsize=16)
    plt.xlabel(ylabel, fontsize=16)
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
    plt.rcParams.update({'font.size': 8})
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.yaxis.get_offset_text().set_fontsize(12) # get the axis offsetText as a Text instance
    plt.text(0.10, 1.1, s, ha='center', va='center', transform=ax.transAxes, fontsize=12)
    plt.text(0.85, 0.5, s0, ha='center', va='center', transform=ax.transAxes, fontsize=14, color='r')
    plt.savefig(savename, dpi = 300, bbox_inches='tight')
    plt.show()

def plot_stl_in_grid_PCA(stl_to_plot, path, param_dict, step_linescan, grid_size, missing_spectra, perverted_spectra, Full, v_min, v_max, Smooth, SVG):
    """Plots STL spectra in a grid using PCA basis.
        
    Args:
        stl_to_plot (str): Filename of the STL data to plot.
        path (list): List containing [datapath2, basispath, basisname, dataname].
        param_dict (dict): Dictionary containing parameters for each spectrum.
        step_linescan (float): Step size for the linescan.
        grid_size (int): Size of the grid (grid_size x grid_size).
        missing_spectra (list): List of indices of missing spectra.
        perverted_spectra (list): List of indices of perverted spectra.
        Full (bool): Whether to plot the full spectrum or a subset.
        v_min (float): Minimum value for the wavelength range.
        v_max (float): Maximum value for the wavelength range.
        Smooth (bool): Whether to apply smoothing to the spectra.
        SVG (bool): If True, save figs in .svg format, else save in .png format
    """

    [datapath2, basispath, basisname, dataname] = path

    date = stl_to_plot[:6]  # Date of acquisition (in the text file)
    spectrum_number = int(stl_to_plot[7:-8])
    print(spectrum_number)

    with open(os.path.join(datapath2, stl_to_plot), 'r', encoding='utf-8') as entete:
        int_time = int(entete.readlines()[5].split()[17][:-1])

    # Load data with cosmic rays already removed

    stl_data = np.load(os.path.join(datapath2, 'Corrected_spectra.npy'))
    wavelength = np.load(os.path.join(datapath2, 'Wavelength.npy'))
    (m, l) = np.shape(stl_data)

    # Get data projected on PCA basis
    coeff = np.load(os.path.join(basispath, dataname))
    vectors = np.load(os.path.join(basispath, basisname))
    for i in range(np.shape(vectors)[0]) :
        vectors[i] = vectors[i]/max(np.abs(vectors[i]))
        if Smooth:
            vectors[i] = savgol_filter(vectors[i], 5, 1)
    
    stl_data = np.matmul(coeff, vectors)
    print(np.shape(stl_data))
    (m, l) = np.shape(stl_data) 

    # Add lines of zeros for missing spectra

    l2 = len(missing_spectra)
    spectra_modif = stl_data[0:m-l2]
    zeroline = np.zeros(l)
    for i in missing_spectra:
        if i > m-l2-1:
            spectra_modif = np.append(spectra_modif, [zeroline], axis=0)
        else:
            spectra_modif = np.insert(spectra_modif, i, zeroline, axis=0)

    # Reshape plots
    if not Full:
        W1 = 1240/v_min
        W2 = 1240/v_max
        i1 = np.argmin(np.abs(wavelength-W1))
        i2 = np.argmin(np.abs(wavelength-W2))

    else:
        i1 = 0
        i2 = -1

    # Plot STL spectra in grid of size: grid_size*grid_size

    for g in range(grid_size):

        num1 = g*grid_size
        num2 = (g+1)*grid_size

        energy = 1240/wavelength
        spectra = spectra_modif[num1:num2]

    # Plot the data: create the 4 figures for plot

        fig2 = plt.figure(2, figsize=[9, 6])
        ax2 = fig2.add_subplot(111, projection='3d')

        fig3 = plt.figure(3, figsize=[9, 6])
        ax3 = fig3.add_subplot(111, projection='3d')

        fig4 = plt.figure(4, figsize=[7, 5])
        ax4 = fig4.add_subplot(111)

        fig5 = plt.figure(5, figsize=[7, 5])
        ax5 = fig5.add_subplot(111)

        for k in range(num2-1, num1-1, -1):

            # remove perverted spectra and lines of zeros added at missing spectra places
            if k in missing_spectra or k in perverted_spectra:
                continue

            j = k - num1
            couleur = j/(num2-num1)
            # cmap = plt.cm.jet
            cmap = plt.cm.rainbow

            # Take one spectrum and remove constant background is needed
            intj = spectra[j][:]

            if Smooth:
                intj = savgol_filter(intj, 5, 1)

        # Create third direction for 3D graph
            Y = np.empty(len(intj))
            Y.fill((j+0.5)*step_linescan)

            # # Corrected data vs wl 3D
            # plt.figure(2)
            # plt.plot(wavelength[i1:i2], Y[i1:i2], Intj[i1:i2],
            #          '-', color=cmap(couleur), linewidth=0.7)

            # Corrected data vs En 3D
            plt.figure(3)
            plt.plot(energy[i1:i2], Y[i1:i2], intj[i1:i2],
                     '-', color=cmap(couleur), linewidth=0.7)

            # # Corrected data vs wl 2D
            # plt.figure(4)
            # plt.plot(wavelength[i1:i2], Intj[i1:i2], '-',
            #          color=cmap(couleur), linewidth=1, label='#{}'.format(Y[0]))

            # Corrected data vs En 2D
            plt.figure(5)
            plt.plot(energy[i1:i2], intj[i1:i2], '-',
                     color=cmap(couleur), linewidth=1, label='#{}'.format(Y[0]))

        # Plotting parameters for graphs and saving
        if spectrum_number in param_dict:
            s = f'line{g} ; {param_dict[spectrum_number][0]}V ; {param_dict[spectrum_number][1]}nA ; VLED = {param_dict[spectrum_number][2]}V ; {int_time}ms'
        else:
            s = f'{int_time}ms'

        foldcontents = os.listdir(basispath)
        foldername = 'Plot_spectra'
        if Smooth:
            foldername += '_Smooth'
        if foldername not in foldcontents:
            os.makedirs(os.path.join(basispath,foldername))
        savepath = os.path.join(basispath,foldername)

        # i = 2
        # ax = ax2
        # xlabel = 'Wavelength (nm)'
        # ylabel = 'Position (nm)'
        # zlabel = 'Intensity (arb. units)'
        # savename = os.path.join(savepath, STL_to_plot.rstrip(
        #     '.txt') + 'vsWl_waterfall' + '_from{}to{}'.format(num1, num2))

        # param_plot_3D(i, ax, xlabel, ylabel, zlabel, s,
        #               date, spectrum_number, savename)

        i = 3
        ax = ax3
        xlabel = 'Energy (eV)'
        ylabel = 'Position (nm)'
        zlabel = 'Intensity (arb. units)'
        if SVG:
            savename = os.path.join(savepath, stl_to_plot.rstrip(
                '.txt') + f'vsEn_waterfall_from{num1}to{num2}.svg')
        else:
            savename = os.path.join(savepath, stl_to_plot.rstrip(
                '.txt') + f'vsEn_waterfall_from{num1}to{num2}')


        param_plot_3d(i, ax, xlabel, ylabel, zlabel, s,
                      date, spectrum_number, savename)

        # i = 4
        # ax = ax4
        # xlabel = 'Wavelength (nm)'
        # ylabel = 'Intensity (arb. units)'
        # savename = os.path.join(savepath, STL_to_plot.rstrip(
        #     '.txt') + 'vsWl_2D' + '_from{}to{}'.format(num1, num2))

        # param_plot_2D(i, ax, xlabel, ylabel, s, date,
        #               spectrum_number, savename)

        i = 5
        ax = ax5
        xlabel = 'Energy (eV)'
        ylabel = 'Intensity (arb. units)'
        if SVG:
            savename = os.path.join(savepath, stl_to_plot.rstrip(
                '.txt') + f'vsEn_2D_from{num1}to{num2}.svg')
        else:
            savename = os.path.join(savepath, stl_to_plot.rstrip(
                '.txt') + f'vsEn_2D_from{num1}to{num2}')

        param_plot_2d(i, ax, xlabel, ylabel, s, date,
                      spectrum_number, savename)
        
        # Save data in text file
        np.savetxt(os.path.join(savepath, stl_to_plot.rstrip('.txt') + f'txt_from{num1}to{num2}.txt'.format(num1, num2)), spectra.T, delimiter='\t')

        plt.close('all')

def plot_sum_spectra(stl_to_plot, datapath2, param_dict, missing_spectra, perverted_spectra, Full, v_min, v_max, Smooth):
    """Plot the sum of the spectra including the FWHM of the main peak

    Args:
        stl_to_plot (str): Text file with the raw acquisition data.
        datapath2 (str): Data path.
        param_dict (dict): Dictoinary with the acquisition parameters.
        missing_spectra (list): List of missing spectra.
        perverted_spectra (list): List of perverted spectra.
        Full (bool): If True plot on the full energy range.
        v_min (float): Minimum value of the energy range.
        v_max (float): Maximum value of the energy range.
        Smooth (bool): If True the data is smoothed.
    """

    spectrum_number = int(stl_to_plot[7:-8])
    print(spectrum_number)

    with open(os.path.join(datapath2, stl_to_plot), 'r', encoding='utf-8') as entete:
        int_time = int(entete.readlines()[5].split()[17][:-1])
    
    # entete = open(os.path.join(datapath2, STL_to_plot), 'r')
    # int_time = int(entete.readlines()[5].split()[17][:-1])

    # Load data with cosmic rays already removed

    stl_data = np.load(os.path.join(datapath2, 'Corrected_spectra.npy'))
    wavelength = np.load(os.path.join(datapath2, 'Wavelength.npy'))
    (m, l) = np.shape(stl_data)

    # Saving folder
    foldcontents = os.listdir(datapath2)
    foldername = 'Plot_spectra'
    if Smooth:
        foldername += '_Smooth'
    if foldername not in foldcontents:
        os.makedirs(os.path.join(datapath2, foldername))
    savepath = os.path.join(datapath2, foldername)

    # Add lines of zeros for missing spectra

    l2 = len(missing_spectra)
    spectra_modif = stl_data[0:m-l2]
    zeroline = np.zeros(l)
    for i in missing_spectra or perverted_spectra:
        if i > m-l2-1:
            spectra_modif = np.append(spectra_modif, [zeroline], axis=0)
        else:
            spectra_modif = np.insert(spectra_modif, i, zeroline, axis=0)

    # Reshape plots
    if not Full:
        W1 = 1240/v_min
        W2 = 1240/v_max
        i1 = np.argmin(np.abs(wavelength-W1))
        i2 = np.argmin(np.abs(wavelength-W2))

    else:
        i1 = 0
        i2 = -1

    energy = 1240/wavelength

    sum_stl = np.sum(stl_data, axis=0)

    ### Finding the FWHM for the summed spectra

    max_spectrum = np.max(sum_stl)
    half_max_spectrum = max_spectrum / 2
    indices = np.where(sum_stl >= half_max_spectrum)[0]
    x_fwhm = energy[indices]
    fwhm = x_fwhm[-1] - x_fwhm[0]

    if spectrum_number in param_dict:
        s = f'{param_dict[spectrum_number][0]}V ; {param_dict[spectrum_number][1]}nA ; VLED = {param_dict[spectrum_number][2]}V ; {int_time}ms'
    else:
        s = f'{int_time}ms'

    plt.plot(energy[i1:i2], sum_stl[i1:i2], 'k', label='summed spectra')
    plt.hlines(half_max_spectrum, x_fwhm[0], x_fwhm[1], color='r', linewidth=2, label='FWHM')
    plt.text((x_fwhm[0] + x_fwhm[-1]) / 2, half_max_spectrum * 1.1, f"{fwhm:.3f}", color="r", fontsize=10, ha="center")

    plt.xlabel('Energy (eV)', fontsize=16)
    plt.ylabel('Intensity (arb. units)', fontsize=16)
    plt.title(s)
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
    plt.rcParams.update({'font.size': 8})

    savename = os.path.join(savepath, stl_to_plot.rstrip('.txt') + '_sum_spectra')
    plt.savefig(savename, dpi=300, bbox_inches='tight')
    
    plt.show()

def plot_sum_spectra_PCA(stl_to_plot, path, param_dict, missing_spectra, perverted_spectra, Full, v_min, v_max, Smooth):
    """Plot the sum of the PCA spectra including the FWHM of the main peak

    Args:
        STL_to_plot (str): data text file
        path (str): filepath where the data is stored
        param_dict (dict): dictionnary with the acquisition parameters
        missing_spectra (list): list of missing spectra
        perverted_spectra (list): list of perverted spectra
        full (bool): plot on the full energy range or not
        Min (int): min energy plotted
        Max (int): max energy plotted
        smooth (bool): data smoothed or not
    """

    [datapath2, basispath, basisname, dataname] = path

    spectrum_number = int(stl_to_plot[7:-8])
    print(spectrum_number)

    with open(os.path.join(datapath2, stl_to_plot), 'r', encoding='utf-8') as entete:
        int_time = int(entete.readlines()[5].split()[17][:-1])

    # Load data with cosmic rays already removed

    stl_data = np.load(os.path.join(datapath2, 'Corrected_spectra.npy'))
    wavelength = np.load(os.path.join(datapath2, 'Wavelength.npy'))
    (m, l) = np.shape(stl_data)

    # Get data projected on PCA basis
    coeff = np.load(os.path.join(basispath, dataname))
    vectors = np.load(os.path.join(basispath, basisname))
    for i in range(np.shape(vectors)[0]) :
        vectors[i] = vectors[i]/max(np.abs(vectors[i]))
        if Smooth:
            vectors[i] = savgol_filter(vectors[i], 5, 1)
    stl_data = np.matmul(coeff, vectors)
    print(np.shape(stl_data))
    (m, l) = np.shape(stl_data)

    # Saving folder

    foldcontents = os.listdir(basispath)
    foldername = 'Plot_spectra'
    if Smooth:
        foldername += '_Smooth'
    if foldername not in foldcontents:
        os.makedirs(os.path.join(basispath,foldername))
    savepath = os.path.join(basispath,foldername)


    # Add lines of zeros for missing spectra

    l2 = len(missing_spectra)
    spectra_modif = stl_data[0:m-l2]
    zeroline = np.zeros(l)
    for i in missing_spectra or perverted_spectra:
        if i > m-l2-1:
            spectra_modif = np.append(spectra_modif, [zeroline], axis=0)
        else:
            spectra_modif = np.insert(spectra_modif, i, zeroline, axis=0)

    # Reshape plots
    if not Full:
        W1 = 1240/v_min
        W2 = 1240/v_max
        i1 = np.argmin(np.abs(wavelength-W1))
        i2 = np.argmin(np.abs(wavelength-W2))

    else:
        i1 = 0
        i2 = -1

    energy = 1240/wavelength

    sum_stl = np.sum(stl_data, axis=0)

    ### Finding the FWHM for the summed spectra
    
    max_spectrum = np.max(sum_stl)
    half_max_spectrum = max_spectrum / 2
    indices = np.where(sum_stl >= half_max_spectrum)[0]
    x_fwhm = energy[indices]
    fwhm = x_fwhm[-1] - x_fwhm[0]

    if spectrum_number in param_dict:
        s = f'{param_dict[spectrum_number][0]}V ; {param_dict[spectrum_number][1]}nA ; VLED = {param_dict[spectrum_number][2]}V ; {int_time}ms'
    else:
        s = f'{int_time}ms'

    plt.plot(energy[i1:i2], sum_stl[i1:i2], 'k', label='summed spectra')
    plt.hlines(half_max_spectrum, x_fwhm[0], x_fwhm[1], color='r', linewidth=2, label='FWHM')
    plt.text((x_fwhm[0] + x_fwhm[-1]) / 2, half_max_spectrum * 1.1, f"{fwhm:.3f}", color="r", fontsize=10, ha="center")

    plt.xlabel('Energy (eV)', fontsize=16)
    plt.ylabel('Intensity (arb. units)', fontsize=16)
    plt.title(s)
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
    plt.rcParams.update({'font.size': 8})

    savename = os.path.join(savepath, stl_to_plot.rstrip('.txt') + '_sum_spectra')
    plt.savefig(savename, dpi=300, bbox_inches='tight')

    plt.xlabel('Energy (eV)', fontsize=16)
    plt.ylabel('Intensity (arb. units)', fontsize=16)
    plt.title(s)
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
    plt.rcParams.update({'font.size': 8})
    plt.legend()

    savename = os.path.join(savepath, stl_to_plot.rstrip('.txt') + '_sum_spectra')
    plt.savefig(savename, dpi=300, bbox_inches='tight')
    
    plt.show()

def integrate_data(DATA_PATH2, PCA, n, Smooth, E1, E2, m, l):
    """Integrate the STLM map spectra and save the result in a txt file.

    Args:
        DATA_PATH2 (str): data path.
        PCA (bool): if True the PCA data are integrated.
        n (int): number of eigen vectors.
        Smooth (bool): if True, the data are smoothed before integration.
        E1 (float): minimum energy for integration.
        E2 (float): maximum energy for integration.
        m (int): number of lines.
        l (int): number of columns.
    """

    ## Use the PCA data
    if PCA:
        name_folder = '{}vecteurs_propres'.format(n)
        basispath = os.path.join(DATA_PATH2, name_folder)
        basisname = "{}vecteurs_propres.npy".format(n)
        dataname = 'coeff_projection_from0to1024.npy'
        # Saving path
        foldcontents = os.listdir(basispath)
        folder_name = 'Correlations'
        if folder_name not in foldcontents:
            os.makedirs(os.path.join(basispath, folder_name))
        saving_path = os.path.join(basispath, folder_name)
        
        ##Get data projected on PCA basis
        coeff = np.load(os.path.join(basispath, dataname))
        vectors = np.load(os.path.join(basispath, basisname))
        for i in range(np.shape(vectors)[0]):
            vectors[i] = vectors[i] / max(vectors[i])
        STL_data = np.matmul(coeff, vectors)

        # (m, l) = np.shape(STL_data)
        # print(m, l)

        ## Get wavelength and energy
        name_wl = 'Wavelength'
        wavelength = np.load(os.path.join(DATA_PATH2, name_wl + '.npy'))
        energy = 1240 / wavelength

    ## Use the raw data
    else:
        ## Get already corrected data
        name_data = 'Corrected_spectra'  # 'Interference_Corrected_spectra' #'Corrected_spectra'
        name_wl = 'Wavelength'
        STL_data = np.load(os.path.join(DATA_PATH2, name_data + '.npy'))
        wavelength = np.load(os.path.join(DATA_PATH2, name_wl + '.npy'))
        energy = 1240 / wavelength
        # Saving path
        foldcontents = os.listdir(DATA_PATH2)
        folder_name = 'Correlations'
        if folder_name not in foldcontents:
            os.makedirs(os.path.join(DATA_PATH2, folder_name))
        saving_path = os.path.join(DATA_PATH2, folder_name)

        # (m, l) = np.shape(STL_data)
        # print(m, l)

        wavelength = np.load(os.path.join(DATA_PATH2, name_wl + '.npy'))
        energy = 1240 / wavelength

    #Bornes intégration
    W1 = 1240 / E1
    W2 = 1240 / E2
    i1 = np.argmin(np.abs(wavelength - W1))
    i2 = np.argmin(np.abs(wavelength - W2))
    #    i1 = 0
    #    i2 = -1

    liste_spectra = np.arange(0, m*l, 1)
    list_intg = []

    for j in liste_spectra:

        print(j)

        spectrum = STL_data[j]

        if Smooth:
            spectrum = savgol_filter(spectrum, 5, 1)

        totalInt = intg.trapz(spectrum[i1:i2], energy[i1:i2])
        list_intg.append(totalInt)

    array_area = np.array(list_intg).reshape(m, l)
    array_area = np.flip(array_area, axis=1)

    np.savetxt(os.path.join(saving_path, f'Intg_matrix_from_{round(1240/W1, 2)}_to_{round(1240/W2, 2)}eV.txt'), array_area, delimiter='\t')