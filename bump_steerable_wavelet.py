import math
import os
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt


######## definition of N, J, L and delta_n
N = 2048 #1024
L = 4 #8
J = 2 #3
max_delta_n = 0
delta_factor = 3;

##### Constant used to compute the filter bank
xi0 = 1.7*np.pi
sigma = 0.248*np.power(2, -0.55)*xi0
c = 1/(1.29) * 2**(L/2 - 1)*math.factorial(L/2-1) / np.sqrt(L/2 * math.factorial(L-2))

# In order to recover Sixin's wavelet, I need to divide hat_phi by this discrepency factor
discrepency_factor = 1.157853617690941



def compute_hat_psi(pulsation_2d, theta, delta_n):
    delta_n_x = delta_n * np.cos(theta)
    delta_n_y = delta_n * np.sin(theta)
    hatpsi = c * np.exp(-np.power(np.abs(pulsation_2d) - xi0, 2) / (xi0**2 - (np.abs(pulsation_2d) - xi0)**2))
    hatpsi = np.nan_to_num(hatpsi)
    hatpsi = hatpsi * (np.absolute(pulsation_2d) <= 2*xi0).astype(int)
    hatpsi = hatpsi * np.power(np.cos(np.angle(pulsation_2d) - theta), L-1)
    hatpsi = hatpsi * (
        (((np.angle(pulsation_2d) - theta) % (2*np.pi)) <= (np.pi / 2.)).astype(int)
        +
        (((np.angle(pulsation_2d) - theta) % (2*np.pi)) >= (3*np.pi / 2.)).astype(int)
    )

    hatpsi = hatpsi*np.exp(- 1j*delta_factor*(np.real(pulsation_2d)*delta_n_x + np.imag(pulsation_2d)*delta_n_y))

    hatpsi[N // 2: N, N // 2: 3 * N // 2] += hatpsi[3 * N // 2:, N // 2: 3 * N // 2]
    hatpsi[N: 3*N//2, N // 2: 3 * N // 2] += hatpsi[:N//2, N // 2: 3 * N // 2]

    hatpsi[N // 2: 3 * N // 2, N // 2: N] += hatpsi[N // 2: 3 * N // 2, 3 * N // 2:]
    hatpsi[N // 2: 3 * N // 2, N: 3 * N // 2] += hatpsi[N // 2: 3 * N // 2, :N // 2]

    hatpsi[N: 3 * N // 2, N:3 * N // 2] += hatpsi[0:N // 2, :N // 2]
    hatpsi[N // 2: N, N // 2:N] += hatpsi[3 * N // 2:, 3 * N // 2:]

    hatpsi[N: 3*N//2, N // 2:N] += hatpsi[:N // 2, 3 * N // 2:]
    hatpsi[N // 2:N, N: 3 * N // 2] += hatpsi[3 * N // 2:, :N // 2]

    return np.fft.ifftshift(hatpsi[N // 2: 3 * N // 2, N // 2: 3 * N // 2])

def compute_hat_phi(puslation_2d):
    return np.exp(- np.absolute(puslation_2d)**2 / (2*(sigma**2)))

def compute_hat_psi_j_theta_delta_n(j, l, delta_n):
    frequencies_x = np.arange(-1, 1, 1/(N)) # np.fft.fftfreq(N)
    frequencies_y = np.arange(-1, 1, 1/(N)) # np.fft.fftfreq(N)
    pulsation_2d = (2**j) * 2 * np.pi * np.array([frequencies_x + 1j * frequency_y for frequency_y in frequencies_y])

    return compute_hat_psi(pulsation_2d, np.pi*l/L, delta_n) # * (2**j)

def compute_hat_phi_j(j):
    frequencies_x = np.fft.fftfreq(N)
    frequencies_y = np.fft.fftfreq(N)
    pulsation_2d = (2**j) * 2 * np.pi * np.array([frequencies_x + 1j * frequency_y for frequency_y in frequencies_y])

    return compute_hat_phi(pulsation_2d) # * (2**j)

def compute_bank_of_wavelet(should_plot = False):
#    matlab_filters = sio.loadmat(
#        './matlab/filters/bumpsteerableg1_fft2d_N' + str(N) + '_J' + str(J) + '_L' + str(L) + '.mat')
#    matlab_fftpsi = matlab_filters['filt_fftpsi']
#    matlab_fftphi = matlab_filters['filt_fftphi']


    list_hat_psi = []
    for l in range(2*L):
        print("l= "+str(l))
        list_hat_psi_l = []
        for j in range(1, J+1):
            list_hat_psi_j = []
            for delta_n in range(max_delta_n+1):
                psi_hat = compute_hat_psi_j_theta_delta_n(j, l, delta_n/2) / discrepency_factor
#                if delta_n == 0 and np.max(np.abs(psi_hat - matlab_fftpsi[j-1,l,:,:])) > 1e-10:
#                    raise Exception('The wavelets created is different from the wavelet from scatnet, while they should be exactly equal.')
                if should_plot and j == 3 and l <= 3:
                    plt.clf()
                    plt.imshow(np.real(np.fft.ifftshift(np.fft.ifft2(psi_hat))))  #, vmin=-0.04, vmax=0.04)
                    plt.colorbar()
                    plt.title("j = " + str(j) + ' l='+str(l) + ' n=' +str(delta_n))
                    plt.show()

                list_hat_psi_j.append(psi_hat)
            list_hat_psi_l.append(np.array(list_hat_psi_j))
        list_hat_psi.append(np.array(list_hat_psi_l))

    hat_phi = compute_hat_phi_j(J-2)
#    if np.max(np.abs(hat_phi-matlab_fftphi)) > 1e-3:  # This matches only at 5e-4, I don't know why.
#        raise Exception

    result = {
        "filt_fftpsi": np.array(list_hat_psi),
        "filt_fftphi": hat_phi
    }
    filename = 'bump_steerable_wavelet_N' + str(N) + '_J' + str(J) + '_L' + str(L) + '_dn' + str(max_delta_n) + '_df' + str(delta_factor) + '.npy'
    np.save(os.path.join('filters', filename), result)

compute_bank_of_wavelet()
