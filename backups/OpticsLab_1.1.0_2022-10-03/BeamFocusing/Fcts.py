#!/usr/bin/env python3.6

import numpy as np
import scipy as sp
import scipy.integrate as integrate
import scipy.special as special
from scipy.optimize import curve_fit
import time as time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import multiprocessing


def fresnel(direction, k1, k2, n1, n2, rho, phi, q):
    kz1 = k1 * np.sqrt(1 - q ** 2)
    kz2 = k2 * np.sqrt(1 - k1 ** 2 / k2 ** 2 * q ** 2)
    ts = 2 * kz1 / (kz1 + kz2)
    tp = 2 * kz1 / (kz1 * (n2 / n1) + kz2 * (n1 / n2))
    x = k1 * rho * q
    if direction == 1:
        fresnel_element = sp.special.jv(0, x) * (ts + tp * kz2 / k2) - \
                          sp.special.jv(2, x) * (tp * kz2 / k2 - ts) * np.cos(2 * phi)
    elif direction == 2:
        fresnel_element = - sp.special.jv(2, x) * (tp * kz2 / k2 - ts) * np.sin(2 * phi)
    elif direction == 3:
        fresnel_element = 1j * sp.special.jv(1, x) * tp * q / k2 * k1 * np.cos(phi)
    else:
        fresnel_element = 0

    return fresnel_element


# By defining an incident field, multiplying it with above fresnel factors and including phase shift in second medium,
# the infinity field E_inf can be calculated. This is then multiplied by additional terms in the so called angular
# spectrum representation (asr).

def asr(direction, k1, k2, n1, n2, f0, num_ap, f, power, z0, rho, phi, z, q):
    kz1 = k1 * np.sqrt(1 - q ** 2)
    kz2 = k2 * np.sqrt(1 - k1 ** 2 / k2 ** 2 * q ** 2)
    # Incident field normalized to a total power of P:
    E_inc = np.sqrt(2 * power / np.pi) / (f0 * num_ap * f) * np.exp(-q ** 2 / (f0 * num_ap) ** 2)
    E_inf = E_inc * np.exp(1j * (kz1 - kz2) * z0) * fresnel(direction, k1, k2, n1, n2, rho, phi, q)
    angular_spectrum_representation = E_inf * (1 - q ** 2) ** (-1 / 4) * k1 / 2 * q * np.exp(
        1j * kz2 * z) * 1j * f * np.exp(-1j * k1 * f)

    return angular_spectrum_representation


# To get the transmitted field, just integrate the asr over the full polar angle covered by the objective. Here a
# substitution is made, such that sin(theta)=q and therefore the bounds are 0 and sin(theta_max)=NA.

def calc_field(direction, k1, k2, n1, n2, f0, num_ap, f, power, z0, rho, phi, z):
    real_integral = integrate.quad(lambda q: np.real(asr(direction, k1, k2, n1, n2, f0, num_ap, f,
                                                         power, z0, rho, phi, z, q)), 0, num_ap)
    imag_integral = integrate.quad(lambda q: np.imag(asr(direction, k1, k2, n1, n2, f0, num_ap, f,
                                                         power, z0, rho, phi, z, q)), 0, num_ap)
    return real_integral[0] + imag_integral[0] * 1j


def field(direction, k1, k2, n1, n2, f0, num_ap, f, power, z0, rho, phi, z):
    transmit_field = np.zeros((np.size(rho), np.size(z)), dtype=np.complex_)
    for iRho, Rho in np.ndenumerate(rho):
        start_time = time.time()
        pool = multiprocessing.Pool(processes=6)
        params = []
        for iZ, Z in np.ndenumerate(z):
            params.append([direction, k1, k2, n1, n2, f0, num_ap, f, power, z0, Rho, phi, Z])
        transmit_field[iRho] = pool.starmap_async(calc_field, params).get()
        pool.close()
        pool.join()
        # transmit_field[iRho, iZ] = calc_field(direction, k1, k2, n1, n2, f0, num_ap, f, power, z0, Rho, phi, Z)
        print('Calculated intensity in {}.direction for z_0 = {:1.0f}nm at radius {}/{}'
              .format(direction, z0 * 10 ** 9, int(iRho[0] + 1), np.size(rho)),
              ' in {:1.3f} sec'.format(time.time() - start_time))
    return transmit_field


# Calculate the focal fields and intensities by the functions defined before:
def calc_intensity(direction, k1, k2, n1, n2, f0, num_ap, f, power, z0, rho, phi, z):
    calc_Field = np.squeeze(field(direction, k1, k2, n1, n2, f0, num_ap, f, power, z0, rho, phi, z))
    calc_Intensity = np.append(np.flipud(np.abs(calc_Field) ** 2), np.abs(calc_Field) ** 2, axis=0) * 10 ** -7
    return calc_Intensity


# Define fct to perform calculation for focus in each spatial point:
def calc_focus(k1, k2, n1, n2, f0, num_ap, f, power, z0, rho, phi, z):
    Intensity = [0, 0, 0]
    for direction in range(1, 4):
        Intensity[direction-1] = calc_intensity(direction, k1, k2, n1, n2, f0, num_ap, f, power, z0, rho, phi, z)
    return Intensity


def calc_focus_axial(wl, k1, k2, n1, n2, f0, num_ap, f, power, z0, rho, phi, z, z_axis):
    Intensity = []
    zpeak = []
    fwhm = []
    for iZ0, Z0 in enumerate(z0):
        Int = np.squeeze(np.abs(field(1, k1, k2, n1, n2, f0, num_ap, f, power, Z0, rho, phi, z))**2) * 10 ** -7
        Intensity.append(Int)
        # Find peak position and FWHM of the intensity via functions defined later:
        zpeak = np.append(zpeak, z_axis[np.where(Int == np.amax(Int))] - Z0 / wl)
        fwhm_axis = z_axis[Int > np.amax(Int) / 2]
        fwhm = np.append(fwhm, np.amax(fwhm_axis) - np.amin(fwhm_axis))
    return {'Intensity': Intensity, 'zPeak': zpeak, 'FWHM': fwhm}


def plot_radial(wl, current_z0, xy_axis, intensity):

    fig2d = plt.figure(figsize=(18, 5))
    fig2d.canvas.set_window_title('Z_0 = {:1.2f} WL'.format(current_z0 / wl))

    for dir_index in range(1, 4):
        axis2d = fig2d.add_subplot(1, 3, dir_index)
        Int = intensity[dir_index - 1][:, np.argmax(intensity[0][int(np.floor(np.size(xy_axis)/2)), :])]
        axis2d.plot(xy_axis, Int)
        plt.tight_layout()
        axis2d.set(xlabel='radial distance / WL', ylabel='intensity $~/~$ kW cm$^{-2}$',
                   ylim=(0, np.amax(Int)))
        dir_label = '|E$_x|^2$' if dir_index == 1 else '|E$_y|^2$' if dir_index == 2 else '|E$_z|^2$'
        axis2d.set_title(dir_label)

    return fig2d


# To plot the intensity vs the xy- and z-axis in a surf-plot, the following function is used. You can provide the
# direction (1,2,3 for x,y,z) of the field, that you want to show as first argument. If you provide 0, all of them will
# be shown.

def plot_2d(wl, z0, xy_axis, z_axis, intensity):
    xyPlane, zPlane = np.meshgrid(xy_axis, z_axis)

    fig2d = plt.figure(figsize=(18, 5))
    for dir_index in range(1, 4):
        axis2d = fig2d.add_subplot(1, 3, dir_index, projection='3d')
        fig2d.canvas.set_window_title("Z_0 = {:1.2f} WL".format(z0 / wl))
        Int = np.squeeze(intensity[dir_index - 1])
        axis2d.plot_surface(xyPlane, zPlane, np.transpose(Int), cmap=cm.jet)
        plt.tight_layout()
        axis2d.set(xlabel='r / WL', ylabel='z / WL', zlabel='intensity $~/~$ kW cm$^{-2}$',
                   zlim=(0, np.amax(Int)))
        dir_label = '|E$_x|^2$' if dir_index == 1 else '|E$_y|^2$' if dir_index == 2 else '|E$_z|^2$'
        axis2d.set_title(dir_label)

    return fig2d


# To plot the focus along z for different z0, you can use:

def plot_axial(z0, z_axis, intensity, n_lines):
    figZ, axis_z = plt.subplots(1, 1, figsize=(18, 5))

    for iz0, Z0 in enumerate(z0):
        # Define a z0 for each line, by picking equaly distant values from z0:
        iz0_shown = np.floor(np.linspace(1, np.size(z0), n_lines)) - 1
        # Plot the intensity vs z-axis curves for each picked z0:
        if iz0 in iz0_shown:
            axis_z.plot(z_axis, intensity[iz0])
            axis_z.set(xlabel='axial position / WL', ylabel='intensity / kW cm$^{-2}$',
                       ylim=(0, np.amax(intensity[iz0])))
            axis_z.set_title('Focus along z for different z$_0$')

    return figZ


# Define linear fit function:
def lin_func(x, a, b):
    return x * a + b


# This function plots and fits you the FWHM for different z0:
def fit_fwhm(z0axis, fwhm):
    fit_params, fit_confs = curve_fit(lin_func, z0axis, fwhm)
    fit_curve = z0axis * fit_params[0] + fit_params[1]
    figFWHM, axisFWHM = plt.subplots(1, 1, figsize=(5, 5))
    axisFWHM.plot(z0axis, fwhm)
    axisFWHM.plot(z0axis, fit_curve)
    axisFWHM.set(xlabel='z$_0$ / WL', ylabel='FWHM / WL', ylim=(0, np.amax(fwhm) * 1.5))
    axisFWHM.set_title('Axial FWHM vs z$_0$')

    return [fit_params, figFWHM]


# This function plots and fits you the peak position for different z0:
def fit_peak(z0axis, z_peak):
    fit_params, fit_confs = curve_fit(lin_func, z0axis, z_peak)
    fit_curve = z0axis * fit_params[0] + fit_params[1]
    fig_peak, axis_peak = plt.subplots(1, 1, figsize=(5, 5))
    axis_peak.plot(z0axis, z_peak)
    axis_peak.plot(z0axis, fit_curve)
    axis_peak.set(xlabel='z$_0$ / WL', ylabel='z$_f$ / WL', ylim=(0, np.amax(z_peak) * 1.5))
    axis_peak.set_title('Axial focus position vs z$_0$')

    return [fit_params, fig_peak]


def evaluate_fwhm(z0_axis, fwhm):
    # Plot and fit the FWHM depending on z0:
    FWHM_fit_params = fit_fwhm(z0_axis, fwhm)[0]
    print('\n   FWHM fit: FWHM = {:1.3f}*z0 + {:1.3f}'.format(float(FWHM_fit_params[0]), float(FWHM_fit_params[1])))
    return FWHM_fit_params


def evaluate_peak_pos(z0_axis, zPeak):
    # Plot and fit the peak position depending on z0:
    PeakPos_fit_params = fit_peak(z0_axis, zPeak)[0]
    print('PeakPos fit: z_f  = {:1.3f}*z0 + {:1.3f}'.format(float(PeakPos_fit_params[0]), float(PeakPos_fit_params[1])))
    return PeakPos_fit_params
