# Imports
import numpy as np
import imageio.v2 as imageio

from datetime import datetime
from functions import *
from parameters import *
import matplotlib.pyplot as plt

# Start the clock
start_time = datetime.now()

# Size of the plot figures
plt.rcParams['figure.figsize'] = [12, 5]

"""
Full routine:
Make use of the Tlusty spectrum and FAMIAS pulsation time-serie
Convolve the 'treated' spectrum with pulsations

Inputs:
Tlusty spectrum
Parameters
FAMIAS pulsation time serie

Output: time serie of spectra affected by pulsations
"""

###############################################
"""Tlusty spectrum scaling and normalisation"""
###############################################

# Import Tlusty spectrum & absolute flux from fils given certain Teff, g, v
Tlusty_spec = Get_Tlusty_spec(common_path, 'Grid Full', Teff, g, v)
Tlusty_sed = Get_Tlusty_spec(common_path, 'Grid SED', Teff, g, v)

# Import the wavelengths for each spectral order
order_numbers, lambda_cen, lambda_min, lambda_max = Get_orders(common_path, order_filename)

# Call theoretical spectrum & SED (scaled with dilution factor)
wavelength, th_flux = Calc_TlustyStellarSpec(Tlusty_spec, rad_sol, dist_pc)
wavelength_SED, th_SED = Calc_TlustyStellarSpec(Tlusty_sed, rad_sol, dist_pc)

# Interpolate the SED over the spectrum wavelength array (don't have the same sampling)
intFunc = interp1d(wavelength_SED, th_SED)
th_SED = intFunc(wavelength)

# Select CubeSPEC wavelength range only
wavelength_cube, th_flux_cube = Cut_specrange(wavelength, th_flux, lambda_lower, lambda_upper)
_, th_SED_cube = Cut_specrange(wavelength, th_SED, lambda_lower, lambda_upper)

# Scale TLusty spectrum to correspond to Vmag given as input
_, th_flux_scaled = Scale_vmag(wavelength_cube, th_flux_cube, Vmag)
_, th_SED_scaled = Scale_vmag(wavelength_cube, th_SED_cube, Vmag)

# Call rebinning function (using the spectrograph )
wave_bin, flux_bin = Rebin(wavelength_cube, th_flux_scaled, stepwidth = mean_wavel/resolution/nPix_per_resol)
_, sed_bin = Rebin(wavelength_cube, th_SED_scaled, stepwidth = mean_wavel/resolution/nPix_per_resol)

# Call broadening mechanisms (to be commented out if done already through FAMIAS)
#flux_rot = Apply_Rotation(wave_bin, flux_bin, vsini) # Rotation
#flux_broad = Apply_Instrument(wave_bin, flux_rot, resolution) # Instrument
#TODO: edge effects on normalised spectrum?! -> seems to be due to instrumental broadening

# Compute the blaze function of the spectrograph
blaze_peak = Calc_Blaze_Peak(wave_bin, lambda_cen, lambda_min, lambda_max, order_numbers)

# Call rebinning & interpolation for necessary parameters (reflectivity, quantum efficiency and transmission)
mirror_reflectivity = Rebin_params(wave_color, mirror_ref_nobin, wave_bin, "quadratic")
QE = Rebin_params(wave_color, QE_nobin, wave_bin, "quadratic")
transmission = Get_transmission(wave_bin, mirror_reflectivity)

# Extract effective flux/SED [photons/A/s]
flux_eff = Calc_ConversionFlux(wave_bin, flux_bin, M1_area, QE, transmission, blaze_peak)
sed_eff = Calc_ConversionFlux(wave_bin, sed_bin, M1_area, QE, transmission, blaze_peak)
# TODO: watch out that the normalisation by the SED gets rid of the blaze_peak shape!
# To keep in mind when returning the observed spectrum (signal + noise)/blaze_peak

# Computes the flux and nominal flux (SED) per spectral bin [photons/bin]
flux_bin, _ = Calc_Flux_bin(wave_bin, flux_eff, ExpTime, resolution, transmission, QE)
sed_bin, _ = Calc_Flux_bin(wave_bin, sed_eff, ExpTime, resolution, transmission, QE)

### Plots the un-normalised spectrum
# plt.figure()
# plt.plot(wave_bin, flux_bin, linewidth=0.7, label='Flux')
# plt.plot(wave_bin, sed_bin, linewidth=0.7, label='SED')
# plt.xlim(lambda_lower - 10, lambda_upper + 10)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Normalised flux')
# plt.legend(loc='best')
# plt.title('Normalised Flux')
# plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini), facecolor='w', transparent=False, format='pdf')
# plt.show()

# Normalise the spectrum
normed_spec = Norm_spec(wave_bin, flux_bin, sed_bin)

# Invert the normalised spectrum
rev_normed_spec = Line_invert(normed_spec)

### Plots the normalised spectrum
# plt.figure()
# plt.plot(wave_bin, normed_spec, linewidth=0.7, label='Normalised Flux')
# plt.plot(wave_bin, rev_normed_spec, linewidth=0.7, label='Inverted Normalised Flux')
# plt.xlim(lambda_lower - 10, lambda_upper + 10)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Normalised flux')
# plt.legend(loc='best')
# plt.title('Normalised Flux')
# plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='pdf')
# plt.show()


###########################################################################
"""Pulsation time serie delivered by FAMIAS: inversion and normalisation"""
###########################################################################

# Import the pulsation profiles
rv, puls = Get_pulsations(pulsation_dir)

# Invert pulsation profiles
rev_puls = Line_invert(puls)

# Compute mean (inverted) pulsation profile
mean_pul_profile = Mean_pul(rev_puls)

# Normalise pulsations (w/ mean or static profile?)
normed_puls = Norm_pul(rv, rev_puls, mean_pul_profile)


#############################################################################
""" Convolution between the normalised spectrum and the pulsation profiles"""
#############################################################################

# Call the convolution between the spectrum and the pulsation profiles (ND array)
convolved_spectra = Convolve_spec_puls(rev_normed_spec, normed_puls, 'same')

### Plot all the convolved spectra
# plt.figure()
# plt.plot(wave_bin, normed_spec + 0.2, linewidth=0.7, label = 'Normalised spectrum')
# for pul in range(np.shape(convolved_spectra)[1]):
#     plt.plot(wave_bin, convolved_spectra[:, pul], linewidth=0.7, label = 'Convolved spectrum')
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Normalised flux')
# plt.axvline(4471.292, lw = 0.5, ls = '--')
# plt.axvline(4861.3, lw = 0.5, ls = '--')
# plt.legend(loc='best')
# plt.show()

#################################################################################
""" Computes the noise & SNR for each spectrum to obtain the Observed spectrum"""
#################################################################################

# New ND array to store the convolved, noisy, normalised spectra to be later written in files
final_time_series = np.zeros_like(convolved_spectra)
# The noise is random, so it needs to be re-computed for each step of the time-series
for step in range(np.shape(final_time_series)[1]):

# However, the noise requires to use the photon flux! So we need to move back to the physical flux again.
    phys_flux_step = convolved_spectra[:, step] * sed_bin

# Compute the straylight for each spectrum of the time-series
    straylight_step = Calc_Stray_only(wave_bin, phys_flux_step, ExpTime, resolution, transmission, QE)

# Use the flux and straylight to compute the noise distribution
    noise_step, snr_step = Calc_Noise(wave_bin, phys_flux_step, straylight_step, ExpTime)

# Combine the actual signal with the noise distribution to obtain the observed signal
    ObsSimFlux = Calc_ObsSimFlux(wave_bin, phys_flux_step, noise_step, blaze_peak)

# Normalise the observed flux
    ObsSimFlux_normed = Norm_spec(wave_bin, ObsSimFlux, sed_bin / blaze_peak)

# Storing the step spectrum in an ND array accounting for the whole time series
    final_time_series[:, step] = ObsSimFlux_normed

########################################################################
""" Write the time series of pulsing spectra as outputs in .dat files"""
########################################################################

# Call the function to save each pulsting spectrum step as .dat files in a given folder
Give_outputs(output_path, 'routine_test_full', wave_bin, final_time_series)

# Stop the clock
end_time = datetime.now()
# Running time
print('Runtime:', end_time - start_time)

#########################################################################
"""Save the plots of the time series and make Gifs animation out of it"""
#########################################################################

filenames = []
for i in range(np.shape(final_time_series)[1]):
    ax = plt.plot(wave_bin, final_time_series[:, i], linewidth = 0.7)
    plt.xlim(4830, 4890)
    plt.ylim(0.58, 1.04)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Normalised flux')
    plt.savefig(common_path+f'/Figures/Gif_test/{i}.png', facecolor='w', transparent=False, dpi = 150)
    filename = common_path+f'/Figures/Gif_test/{i}.png'
    filenames.append(filename)
    plt.close()

filenames_pul = []
for i in range(np.shape(puls)[1]):
    ax = plt.plot(rv, puls[:, i], linewidth = 0.7)
    plt.xlabel('RV [km/s]')
    plt.ylabel('Normalised flux')
    plt.ylim(0.94, 1.01 )
    plt.savefig(common_path+f'/Figures/Gif_pul_test/{i}.png', facecolor='w', transparent=False, dpi = 150)
    filename_pul = common_path+f'/Figures/Gif_pul_test/{i}.png'
    filenames_pul.append(filename_pul)
    plt.close()

filenames_mix = []
for i in range(np.shape(final_time_series)[1]):
    fig, ax = plt.subplots(2, 1, dpi = 150, figsize = (8, 10))
    ax[0].plot(wave_bin, final_time_series[:, i], linewidth=0.7)
    ax[0].set_xlim(4830, 4890)
    ax[0].set_ylim(0.58, 1.04)
    ax[0].set_xlabel('Wavelength [$\AA$]')
    ax[0].set_ylabel('Normalised flux')
    ax[0].set_title('LPV')
    #
    ax[1].plot(rv, puls[:, i], linewidth=0.7)
    ax[1].set_xlabel('RV [km/s]')
    ax[1].set_ylabel('Normalised flux')
    ax[1].set_ylim(0.94, 1.01)
    ax[1].set_title('Pulsation profiles')
    #
    plt.savefig(common_path + f'/Figures/Gif_mix_test/{i}.png', facecolor='w', transparent=False, dpi=150)
    filename_mix = common_path + f'/Figures/Gif_mix_test/{i}.png'
    filenames_mix.append(filename_mix)
    # plt.tight_layout()
    plt.close()

### Creating the Gifs from previously saved figures
with imageio.get_writer(common_path+'/Figures/myAnimation_vsini100.gif', fps=20) as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
print('Line Gif saved!')

with imageio.get_writer(common_path+'/Figures/myAnimation_pul.gif', fps=20) as writer:
    for filename in filenames_pul:
        image = imageio.imread(filename)
        writer.append_data(image)
print('Pul Gif saved!')

with imageio.get_writer(common_path+'/Figures/myAnimation_mix.gif', fps=20) as writer:
    for filename in filenames_mix:
        image = imageio.imread(filename)
        writer.append_data(image)
print('Mix Gif saved!')