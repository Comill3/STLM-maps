import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter, find_peaks, peak_widths


def param_plot_2D(i, ax, xlabel, ylabel, s, date, spectrum_number, savename):
    """Plot parameters for 2D spectra

    Args:
        i (int): _description_
        ax (Axes): fig.add_subplot()
        xlabel (string): x axis name
        ylabel (string): y axis name
        s (string): acquisition parameters
        date (string): _description_
        spectrum_number (int): Acquisition number
        savename (string): fig_name.png
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


def param_plot_3D(i, ax, xlabel, ylabel, zlabel, s, date, spectrum_number, savename):
    """Plot parameters for 3D spectra

    Args:
        i (int): _description_
        ax (Axes): fig.add_subplot()
        xlabel (string): x axis name
        ylabel (string): y axis name
        zlabel (string): z axis name
        s (string): acquisition parameters
        date (string): _description_
        spectrum_number (int): Acquisition number
        savename (string): fig_name.png
    """
    plt.figure(i)
    ax.set_zlabel(zlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=12)
    # Change the backgroung color
    ax.xaxis.set_pane_color((0.5, 0.5, 0.5, 1.0))
    ax.yaxis.set_pane_color((0.5, 0.5, 0.5, 1.0))
    ax.zaxis.set_pane_color((0.5, 0.5, 0.5, 1.0))

    # Additional text on the spectra
    # place text on a fixed position on the ax object
    ax.text2D(0.4, 0.95, s, transform=ax.transAxes, fontsize=10)
    ax.text2D(0.15, 0.95, str(date) + ' #' + str(spectrum_number),
              transform=ax.transAxes, fontsize=10)

    plt.rcParams.update({'font.size': 8})
    # Set the elevation and azimuth of the axes in degrees (not radians)
    ax.view_init(35, -80)
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=10)
    ax.tick_params(axis='z', labelsize=10)
    plt.savefig(savename, dpi=300, bbox_inches='tight')


def plot_STL_in_grid(STL_to_plot, datapath2, param_dict, step_linescan, grid_size, missing_spectra, perverted_spectra, full, min, max, smooth):
    """Plotting the spectra from the same line together in 3D and in 2D, also it saves them in a specific folder

    Args:
        STL_to_plot (string): acquisition text file
        datapath2 (string): path where the acquisition text file is
        param_dict (dictionnary): _dictionnary with acquisition parameters
        step_linescan (float): line size/number of spectra per line
        grid_size (int): number of spectra per line
        missing_spectra (list): list of missing spectra
        perverted_spectra (list): list of missing spectra + 1
        full (boolean): ask if we plot on the full range or not
        min (float): min energy plotted
        max (float): max energy plotted
        smooth (boolean): ask if the spectra has to be smoothed or not
    """

    date = STL_to_plot[:6]  # Date of acquisition (in the text file)
    spectrum_number = int(STL_to_plot[7:-8])
    print(spectrum_number)

    entete = open(os.path.join(datapath2, STL_to_plot), 'r')
    int_time = int(entete.readlines()[5].split()[17][:-1])

    # Load data with cosmic rays already removed

    name_data = 'Corrected_spectra' #'Corrected_spectra' #'Interference_Corrected_spectra'
    name_wl = 'Wavelength'

    STL_data = np.load(os.path.join(datapath2, name_data + '.npy'))
    wavelength = np.load(os.path.join(datapath2, name_wl + '.npy'))
    (m, l) = np.shape(STL_data)

    # Add lines of zeros for missing spectra

    l2 = len(missing_spectra)
    Spectra_modif = STL_data[0:m-l2]
    zeroline = np.zeros(l)
    for i in missing_spectra:
        if i > m-l2-1:
            Spectra_modif = np.append(Spectra_modif, [zeroline], axis=0)
        else:
            Spectra_modif = np.insert(Spectra_modif, i, zeroline, axis=0)

    # Reshape plots
    if full == False:
        W1 = 1240/min
        W2 = 1240/max
        i1 = np.argmin(np.abs(wavelength-W1))
        i2 = np.argmin(np.abs(wavelength-W2))

    elif full == True:
        i1 = 0
        i2 = -1

    # Plot STL spectra in grid of size: grid_size*grid_size

    for g in range(grid_size):

        num1 = g*grid_size
        num2 = (g+1)*grid_size

        energy = 1240/wavelength
        Spectra = Spectra_modif[num1:num2]

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
            Intj = Spectra[j][:]

            if smooth:
                Intj = savgol_filter(Intj, 5, 1)

        # Create third direction for 3D graph
            Y = np.empty(len(Intj))
            Y.fill((j+0.5)*step_linescan)

            # Corrected data vs wl 3D
            plt.figure(2)
            plt.plot(wavelength[i1:i2], Y[i1:i2], Intj[i1:i2],
                     '-', color=cmap(couleur), linewidth=0.7)

            # Corrected data vs En 3D
            plt.figure(3)
            plt.plot(energy[i1:i2], Y[i1:i2], Intj[i1:i2],
                     '-', color=cmap(couleur), linewidth=0.7)

            # Corrected data vs wl 2D
            plt.figure(4)
            plt.plot(wavelength[i1:i2], Intj[i1:i2], '-',
                     color=cmap(couleur), linewidth=1, label='#{}'.format(Y[0]))

            # Corrected data vs En 2D
            plt.figure(5)
            plt.plot(energy[i1:i2], Intj[i1:i2], '-',
                     color=cmap(couleur), linewidth=1, label='#{}'.format(Y[0]))

        # Plotting parameters for graphs and saving
        if spectrum_number in param_dict:
            s = 'line{} ; {}V ; {}nA ; VLED = {}V ; {}ms'.format(g, param_dict[spectrum_number][0], param_dict[spectrum_number][1], param_dict[spectrum_number][2], int_time)
        else:
            s = '{}ms'.format(int_time)

        foldcontents = os.listdir(datapath2)
        foldername = 'Plot_spectra'
        if smooth == True:
            foldername += '_Smooth'
        if foldername not in foldcontents:
            os.makedirs(os.path.join(datapath2, foldername))
        savepath = os.path.join(datapath2, foldername)

        i = 2
        ax = ax2
        xlabel = 'Wavelength (nm)'
        ylabel = 'Position (nm)'
        zlabel = 'Intensity (arb. units)'
        savename = os.path.join(savepath, STL_to_plot.rstrip(
            '.txt') + 'vsWl_waterfall' + '_from{}to{}'.format(num1, num2))

        param_plot_3D(i, ax, xlabel, ylabel, zlabel, s,
                      date, spectrum_number, savename)

        i = 3
        ax = ax3
        xlabel = 'Energy (eV)'
        ylabel = 'Position (nm)'
        zlabel = 'Intensity (arb. units)'
        savename = os.path.join(savepath, STL_to_plot.rstrip(
            '.txt') + 'vsEn_waterfall' + '_from{}to{}'.format(num1, num2))

        param_plot_3D(i, ax, xlabel, ylabel, zlabel, s,
                      date, spectrum_number, savename)

        i = 4
        ax = ax4
        xlabel = 'Wavelength (nm)'
        ylabel = 'Intensity (arb. units)'
        savename = os.path.join(savepath, STL_to_plot.rstrip(
            '.txt') + 'vsWl_2D' + '_from{}to{}'.format(num1, num2))

        param_plot_2D(i, ax, xlabel, ylabel, s, date,
                      spectrum_number, savename)

        i = 5
        ax = ax5
        xlabel = 'Energy (eV)'
        ylabel = 'Intensity (arb. units)'
        savename = os.path.join(savepath, STL_to_plot.rstrip(
            '.txt') + 'vsEn_2D' + '_from{}to{}'.format(num1, num2))

        param_plot_2D(i, ax, xlabel, ylabel, s, date,
                      spectrum_number, savename)

        plt.close('all')


def plot_one_spectrum(X, Y, couleur, xlabel, ylabel, s, s0, savename, ax):
        """Plot parameters for one unique spectra

        Args:
            X (array): X data
            Y (array): Y data
            couleur (string): color of the plot
            xlabel (string): x axis name
            ylabel (string): y axis name
            s (string): text
            s0 (string): text
            savename (string): savename.png
            ax (Axes): fig.add_subplot()
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



def plot_STL_in_grid_PCA(STL_to_plot, path, param_dict, step_linescan, grid_size, missing_spectra, perverted_spectra, full, Min, Max, smooth):
    """Plotting the spectra from the same line together in 3D and in 2D, also it saves them in a specific folder

    Args:
        STL_to_plot (string): acquisition text file
        path (string): path where the acquisition text file is
        param_dict (dictionnary): _dictionnary with acquisition parameters
        step_linescan (float): line size/number of spectra per line
        grid_size (int): number of spectra per line
        missing_spectra (list): list of missing spectra
        perverted_spectra (list): list of missing spectra + 1
        full (boolean): ask if we plot on the full range or not
        Min (float): min energy plotted
        Max (float): max energy plotted
        smooth (boolean): ask if the spectra has to be smoothed or not
    """

    [datapath2, basispath, basisname, dataname] = path

    date = STL_to_plot[:6]  # Date of acquisition (in the text file)
    spectrum_number = int(STL_to_plot[7:-8])
    print(spectrum_number)

    entete = open(os.path.join(datapath2, STL_to_plot), 'r')
    int_time = int(entete.readlines()[5].split()[17][:-1])

    # Load data with cosmic rays already removed

    STL_data = np.load(os.path.join(datapath2, 'Corrected_spectra.npy'))
    wavelength = np.load(os.path.join(datapath2, 'Wavelength.npy'))
    (m, l) = np.shape(STL_data)

    # Get data projected on PCA basis
    coeff = np.load(os.path.join(basispath, dataname))
    vectors = np.load(os.path.join(basispath, basisname))
    for i in range(np.shape(vectors)[0]) :
        vectors[i] = vectors[i]/max(np.abs(vectors[i]))
        if smooth == True :
            vectors[i] = savgol_filter(vectors[i], 5, 1)
    STL_data = np.matmul(coeff, vectors)
    print(np.shape(STL_data))
    (m, l) = np.shape(STL_data) 

    # Add lines of zeros for missing spectra

    l2 = len(missing_spectra)
    Spectra_modif = STL_data[0:m-l2]
    zeroline = np.zeros(l)
    for i in missing_spectra:
        if i > m-l2-1:
            Spectra_modif = np.append(Spectra_modif, [zeroline], axis=0)
        else:
            Spectra_modif = np.insert(Spectra_modif, i, zeroline, axis=0)

    # Reshape plots
    if full == False:
        W1 = 1240/Min
        W2 = 1240/Max
        i1 = np.argmin(np.abs(wavelength-W1))
        i2 = np.argmin(np.abs(wavelength-W2))

    elif full == True:
        i1 = 0
        i2 = -1

    # Plot STL spectra in grid of size: grid_size*grid_size

    for g in range(grid_size):

        num1 = g*grid_size
        num2 = (g+1)*grid_size

        Energy = 1240/wavelength
        Spectra = Spectra_modif[num1:num2]

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
            Intj = Spectra[j][:]

            if smooth:
                Intj = savgol_filter(Intj, 5, 1)

        # Create third direction for 3D graph
            Y = np.empty(len(Intj))
            Y.fill((j+0.5)*step_linescan)

            # # Corrected data vs wl 3D
            # plt.figure(2)
            # plt.plot(wavelength[i1:i2], Y[i1:i2], Intj[i1:i2],
            #          '-', color=cmap(couleur), linewidth=0.7)

            # Corrected data vs En 3D
            plt.figure(3)
            plt.plot(Energy[i1:i2], Y[i1:i2], Intj[i1:i2],
                     '-', color=cmap(couleur), linewidth=0.7)

            # # Corrected data vs wl 2D
            # plt.figure(4)
            # plt.plot(wavelength[i1:i2], Intj[i1:i2], '-',
            #          color=cmap(couleur), linewidth=1, label='#{}'.format(Y[0]))

            # Corrected data vs En 2D
            plt.figure(5)
            plt.plot(Energy[i1:i2], Intj[i1:i2], '-',
                     color=cmap(couleur), linewidth=1, label='#{}'.format(Y[0]))

        # Plotting parameters for graphs and saving
        if spectrum_number in param_dict:
            s = 'line{} ; {}V ; {}nA ; VLED = {}V ; {}ms'.format(g, param_dict[spectrum_number][0], param_dict[spectrum_number][1], param_dict[spectrum_number][2], int_time)
        else:
            s = '{}ms'.format(int_time)

        foldcontents = os.listdir(basispath)
        foldername = 'Plot_spectra'
        if smooth == True :
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
        savename = os.path.join(savepath, STL_to_plot.rstrip(
            '.txt') + 'vsEn_waterfall' + '_from{}to{}'.format(num1, num2))

        param_plot_3D(i, ax, xlabel, ylabel, zlabel, s,
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
        savename = os.path.join(savepath, STL_to_plot.rstrip(
            '.txt') + 'vsEn_2D' + '_from{}to{}'.format(num1, num2))

        param_plot_2D(i, ax, xlabel, ylabel, s, date,
                      spectrum_number, savename)
        
        # Save data in text file
        np.savetxt(os.path.join(savepath, STL_to_plot.rstrip('.txt') + 'txt' + '_from{}to{}.txt'.format(num1, num2)), Spectra.T, delimiter='\t')

        plt.close('all')

def plot_sum_spectra(STL_to_plot, datapath2, param_dict, missing_spectra, perverted_spectra, full, Min, Max, smooth):
    """Plot the sum of the raw spectra including the FWHM of the main peak

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

    spectrum_number = int(STL_to_plot[7:-8])
    print(spectrum_number)

    entete = open(os.path.join(datapath2, STL_to_plot), 'r')
    int_time = int(entete.readlines()[5].split()[17][:-1])

    # Load data with cosmic rays already removed

    STL_data = np.load(os.path.join(datapath2, 'Corrected_spectra.npy'))
    wavelength = np.load(os.path.join(datapath2, 'Wavelength.npy'))
    (m, l) = np.shape(STL_data)

    # Saving folder

    foldcontents = os.listdir(datapath2)
    foldername = 'Plot_spectra'
    if smooth == True:
        foldername += '_Smooth'
    if foldername not in foldcontents:
        os.makedirs(os.path.join(datapath2, foldername))
    savepath = os.path.join(datapath2, foldername)



    # Add lines of zeros for missing spectra

    l2 = len(missing_spectra)
    Spectra_modif = STL_data[0:m-l2]
    zeroline = np.zeros(l)
    for i in missing_spectra or perverted_spectra:
        if i > m-l2-1:
            Spectra_modif = np.append(Spectra_modif, [zeroline], axis=0)
        else:
            Spectra_modif = np.insert(Spectra_modif, i, zeroline, axis=0)

    # Reshape plots
    if full == False:
        W1 = 1240/Min
        W2 = 1240/Max
        i1 = np.argmin(np.abs(wavelength-W1))
        i2 = np.argmin(np.abs(wavelength-W2))

    elif full == True:
        i1 = 0
        i2 = -1

    energy = 1240/wavelength

    sum_STL = np.sum(STL_data, axis=0)

    ### Removing the missing spectra


    ### Finding the FWHM for the summed spectra
    
    # Finding the peaks
    peaks, _ = find_peaks(sum_STL, height=500000, distance=100000)
    # Calculate the fwhm
    results_half = peak_widths(sum_STL, peaks, rel_height=0.5)
    print(results_half)

    fwhm_values = results_half[0]
    fwhm_positions = results_half[2:4]

    if spectrum_number in param_dict:
        s = '{}V ; {}nA ; VLED = {}V ; {}ms'.format(param_dict[spectrum_number][0], param_dict[spectrum_number][1], param_dict[spectrum_number][2], int_time)
    else:
        s = '{}ms'.format(int_time)

    plt.plot(energy[i1:i2], sum_STL[i1:i2], 'k', label='summed spectra')
    plt.plot(energy[peaks], sum_STL[peaks], 'x')

    for peak, fwhm, (left_ips, right_ips) in zip(peaks, fwhm_values, zip(*fwhm_positions)):
        left_base = energy[int(left_ips)]
        right_base = energy[int(right_ips)]
        plt.hlines(sum_STL[peak] / 2, left_base, right_base, color='r', label='FWHM' if peak == peaks[0] else "")
        plt.text(right_base, sum_STL[peak]/2, f'{       (right_base - left_base)*10**3:.2f} meV', color='r', horizontalalignment='left', verticalalignment='bottom')

    plt.xlabel('Energy (eV)', fontsize=16)
    plt.ylabel('Intensity (arb. units)', fontsize=16)
    plt.title(s)
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
    plt.rcParams.update({'font.size': 8})
    plt.legend()

    savename = os.path.join(savepath, STL_to_plot.rstrip('.txt') + '_sum_spectra')
    plt.savefig(savename, dpi=300, bbox_inches='tight')
    
    plt.show()

def plot_sum_spectra_PCA(STL_to_plot, path, param_dict, missing_spectra, perverted_spectra, full, Min, Max, smooth):
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

    spectrum_number = int(STL_to_plot[7:-8])
    print(spectrum_number)

    entete = open(os.path.join(datapath2, STL_to_plot), 'r')
    int_time = int(entete.readlines()[5].split()[17][:-1])

    # Load data with cosmic rays already removed

    STL_data = np.load(os.path.join(datapath2, 'Corrected_spectra.npy'))
    wavelength = np.load(os.path.join(datapath2, 'Wavelength.npy'))
    (m, l) = np.shape(STL_data)

    # Get data projected on PCA basis
    coeff = np.load(os.path.join(basispath, dataname))
    vectors = np.load(os.path.join(basispath, basisname))
    for i in range(np.shape(vectors)[0]) :
        vectors[i] = vectors[i]/max(np.abs(vectors[i]))
        if smooth == True :
            vectors[i] = savgol_filter(vectors[i], 5, 1)
    STL_data = np.matmul(coeff, vectors)
    print(np.shape(STL_data))
    (m, l) = np.shape(STL_data)

    # Saving folder

    foldcontents = os.listdir(basispath)
    foldername = 'Plot_spectra'
    if smooth == True :
        foldername += '_Smooth'
    if foldername not in foldcontents:
        os.makedirs(os.path.join(basispath,foldername))
    savepath = os.path.join(basispath,foldername)


    # Add lines of zeros for missing spectra

    l2 = len(missing_spectra)
    Spectra_modif = STL_data[0:m-l2]
    zeroline = np.zeros(l)
    for i in missing_spectra or perverted_spectra:
        if i > m-l2-1:
            Spectra_modif = np.append(Spectra_modif, [zeroline], axis=0)
        else:
            Spectra_modif = np.insert(Spectra_modif, i, zeroline, axis=0)

    # Reshape plots
    if full == False:
        W1 = 1240/Min
        W2 = 1240/Max
        i1 = np.argmin(np.abs(wavelength-W1))
        i2 = np.argmin(np.abs(wavelength-W2))

    elif full == True:
        i1 = 0
        i2 = -1

    energy = 1240/wavelength

    sum_STL = np.sum(STL_data, axis=0)

    ### Removing the missing spectra


    ### Finding the FWHM for the summed spectra
    
    # Finding the peaks
    peaks, _ = find_peaks(sum_STL, height=500000, distance=100000)
    # Calculate the fwhm
    results_half = peak_widths(sum_STL, peaks, rel_height=0.5)
    print(results_half)

    fwhm_values = results_half[0]
    fwhm_positions = results_half[2:4]

    if spectrum_number in param_dict:
        s = '{}V ; {}nA ; VLED = {}V ; {}ms'.format(param_dict[spectrum_number][0], param_dict[spectrum_number][1], param_dict[spectrum_number][2], int_time)
    else:
        s = '{}ms'.format(int_time)

    plt.plot(energy[i1:i2], sum_STL[i1:i2], 'k', label='summed spectra')
    plt.plot(energy[peaks], sum_STL[peaks], 'x')

    for peak, fwhm, (left_ips, right_ips) in zip(peaks, fwhm_values, zip(*fwhm_positions)):
        left_base = energy[int(left_ips)]
        right_base = energy[int(right_ips)]
        plt.hlines(sum_STL[peak] / 2, left_base, right_base, color='r', label='FWHM' if peak == peaks[0] else "")
        plt.text(right_base, sum_STL[peak]/2, f'{       (right_base - left_base)*10**3:.2f} meV', color='r', horizontalalignment='left', verticalalignment='bottom')

    plt.xlabel('Energy (eV)', fontsize=16)
    plt.ylabel('Intensity (arb. units)', fontsize=16)
    plt.title(s)
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0))
    plt.rcParams.update({'font.size': 8})
    plt.legend()

    savename = os.path.join(savepath, STL_to_plot.rstrip('.txt') + '_sum_spectra')
    plt.savefig(savename, dpi=300, bbox_inches='tight')
    
    plt.show()