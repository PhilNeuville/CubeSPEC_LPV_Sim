# Imports
import os
from os.path import isfile
import sys
import six
import json

import pandas as pd
import numpy as np
import inspect
import matplotlib.pyplot as plt
from scipy import signal

# Import scripts
from functions import *
from parameters import *

# Import Tlusty spectrum from file given certain Teff, g, v
Tlusty_spec = get_Tlusty_spec(common_path, 'Grid Full', Teff, g, v)
Tlusty_sed = get_Tlusty_spec(common_path, 'Grid SED', Teff, g, v)

# Import the wavelengths for each order
order_numbers, lambda_cen, lambda_min, lambda_max = get_orders(common_path, order_filename)

# Settings for figures
#plt.style.use('./Plots/plots.mplstyle')
# Check whether the specified path exists or not
isExist = os.path.exists(fig_path)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(fig_path)
    print("The new directory is created!")

"""
PART 1
Code to caculate observed simulated spectrum O(lambda) = t(lambda) * r(lambda) + n(lambda)
with:
t(lambda) = Theoretical flux: flux (from Tlusty) * R²/D²;
r(lambda) = Response curve;
n(lambda) = Noise; Gaussian/Normal noise: m=0, sigma=sqrt(F)
"""

# Call theoretical spectrum (scaled with dilution factor)
wavelength, th_flux = Calc_TlustyStellarSpec(Tlusty_spec, rad_sol, dist_pc)
wavelength_SED, th_SED = Calc_TlustyStellarSpec(Tlusty_sed, rad_sol, dist_pc)

# Interpolate the SED over the spectrum wavelength range (don't have the same x-axis)
intFunc = interp1d(wavelength_SED, th_SED)
th_SED = intFunc(wavelength)

print(len(wavelength), len(th_SED))

# Select only CubeSPEC wavelength range
wavelength_cube, th_flux_cube = cut_specrange(wavelength, th_flux, lambda_lower, lambda_upper)
_, th_SED_cube = cut_specrange(wavelength, th_SED, lambda_lower, lambda_upper)

# Size of the plot figures
plt.rcParams['figure.figsize'] = [12, 5]

# Plot theoretical spectrum + highlight CubeSPEC range
# plt.figure()
# fig, ax = plt.subplots()
# ax.plot(wavelength, th_flux, linewidth=0.4, label='Full range')
# ax.scatter(wavelength_SED, th_SED, color = 'orange', label='Full range SED', s = 5)
# ax.plot(wavelength_cube, th_flux_cube, linewidth=0.5, color='red', label='CubeSPEC range')
# ax.scatter(wavelength_SED_cube, th_SED_cube, color = 'green', label='Full range SED', s = 5)
# axin = ax.inset_axes([0.48, 0.4, 0.5, 0.4])
# axin.plot(wavelength_cube, th_flux_cube, linewidth=0.5, color='red')
# axin.set_xlim(lambda_lower, lambda_upper)
# ax.indicate_inset_zoom(axin)
# ax.set_xlabel('Wavelength [$\AA$]')
# ax.set_ylabel('Flux [erg/s/cm²/$\AA$]')
# ax.legend(loc='best')
# ax.set_title('Tlusty theoretical stellar spectrum')
# #plt.savefig(fig_path + 'TlustyTheoSpectrum_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini), facecolor='w', transparent=False, format='svg')
# plt.show()

# Plot CubeSPEC theoretical spectrum range only
# plt.figure()
# plt.plot(wavelength_cube, th_flux_cube, linewidth=0.5, label='CubeSpec range')
# plt.xlim(lambda_lower, lambda_upper)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Flux [erg/s/cm²/$\AA$]')
# plt.legend(loc='best')
# plt.title('Tlusty stellar spectrum (CubeSpec range)')
# plt.savefig(fig_path + 'TheoreticalSpectrum_T{}_D{}_R{}_g{}.svg'.format(Teff, dist_pc, rad_sol, g),
#             facecolor='w', transparent=False, format='svg')
# # plt.show()

# Estimate V-band magnitude
tcs = pyasl.TransmissionCurves()
tcs.addSpitzerIRACPassbands()
print("Available bands: ", tcs.availableBands())

# Plot transmission curves for Johnson U, B, and V bands
# plt.figure()
# for (b, c) in six.iteritems({"U": "m", "B": "b", "V": "k"}):
#     tc = tcs.getTransCurve("Johnson " + b)
#     trans = tc(wavelength_cube)
#     plt.plot(wavelength_cube, trans, c + '--', label="Johnson " + b)
#
# plt.title('Transmission Curves for Johnson U, B, V bands')
# plt.legend()
# plt.xlim(lambda_lower, lambda_upper)
# plt.xlabel("Wavelength [$\AA$]")
# plt.ylabel("Transmission")
# plt.show()

# Multiply spectrum by transmission curve for V band = Filter out all the flux that falls outside the V band
tc = tcs.getTransCurve("Johnson V")
trans = tc(wavelength_cube)

# Integrate product of spectrum and transmission curve over all wavelengths = Give total flux within the V band
flux_in_v_band = simps(th_flux_cube * trans, wavelength_cube)
flux_in_v_band_SED = simps(th_SED_cube * trans, wavelength_cube)
# flux_in_v_band = simps(th_flux_cube[indices], wavelength_cube[indices])

# Convert flux in V band to magnitude using zero-point flux for V band
v_mag = -2.5 * np.log10(flux_in_v_band / flux_0_johnson_v)
v_mag_SED = -2.5 * np.log10(flux_in_v_band_SED / flux_0_johnson_v)

print("Estimated V-band magnitude before scaling: ", v_mag, v_mag_SED)

# Scale Tlusty spectrum to correspond to V magnitude given as input
_, th_flux_scaled = scale_vmag(wavelength_cube, th_flux_cube, Vmag)
_, th_SED_scaled = scale_vmag(wavelength_cube, th_SED_cube, Vmag)

# Estimation of V-band magnitude after scaling (Check)
# Integrate product of spectrum and transmission curve over all wavelengths = Give total flux within the V band
flux_in_v = simps(th_flux_scaled * trans, wavelength_cube)
flux_in_v_SED = simps(th_SED_scaled * trans, wavelength_cube)
# flux_in_v = simps(th_flux_scaled[indices], wavelength_cube[indices])

# Convert flux in V band to magnitude using zero-point flux for V band
v_mag_sc = -2.5 * np.log10(flux_in_v / flux_0_johnson_v)
v_mag_sc_SED = -2.5 * np.log10(flux_in_v_SED / flux_0_johnson_v)

print("Estimated V-band magnitude after scaling: ", v_mag_sc, v_mag_sc_SED)

# Plot CubeSPEC scaled flux
# plt.figure()
# plt.plot(wavelength_cube, th_flux_scaled, linewidth=0.4, label='Flux values')
# plt.xlim(lambda_lower, lambda_upper)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Flux [erg/s/cm²/$\AA$]')
# plt.title('Scaled theoretical flux to Vmag = {}'.format(Vmag))
# plt.savefig(fig_path + 'ScaledSpectrum_Vmag{}_T{}g{}D{}R{}vsini{}.svg'.format(Vmag, Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='svg')
# plt.legend()
# plt.show()

#mean_wavel = 5000        # middle wavelength of CubeSPEC range
#resol_power = resolution # not the same as the resolution dlambda!!! R = lambda / dlambda = 55000 is the resolving power

# Call rebinning function
wave_bin, flux_bin = rebin(wavelength_cube, th_flux_scaled, stepwidth = mean_wavel/resolution/nPix_per_resol)
_, sed_bin = rebin(wavelength_cube, th_SED_scaled, stepwidth = mean_wavel/resolution/nPix_per_resol)

# Call rotational broadening mechanism
flux_rot = Apply_Rotation(wave_bin, flux_bin, vsini)
flux_broad = Apply_Instrument(wave_bin, flux_rot, resolution)

# Plot line-broadened theoretical spectrum
# plt.figure()
# plt.plot(wave_bin, flux_broad, linewidth = 0.7, label='Flux values')
# plt.xlim(lambda_lower, lambda_upper)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Flux [erg/s/cm²/$\AA$]')
# plt.legend(loc='best')
# plt.title('Line broadened theoretical spectrum')
# plt.savefig(fig_path+'LineBroadenedTheoSpectrum_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='svg')
# plt.show()

# Compare line-broadened & not line-broadened spectrum
# plt.figure()
# plt.plot(wave_bin, flux_broad, linewidth = 0.8, label='Flux broad',color='red')
# plt.scatter(wave_bin, sed_bin, color = 'orange', label = 'Absolute flux', s = 5)
# plt.plot(wavelength_cube, th_flux_scaled, linewidth = 0.5, alpha=0.5, label='Flux no broad')
# plt.xlim(lambda_lower, lambda_upper)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Flux [erg/s/cm²/$\AA$]')
# plt.legend(loc='best')
# plt.title('Theoretical vs. Line broadened')
# plt.savefig(fig_path+'TheoLineBroadSpec_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='svg')
# plt.show()

line_lower = 4840
line_upper = 4880

indices2 = np.where((wave_bin >= line_lower) & (wave_bin < line_upper))
indices3 = np.where((wavelength_cube >= line_lower) & (wavelength_cube < line_upper))

# Compare line-broadened & not line-broadened spectrum - ZOOM
# plt.figure()
# plt.plot(wave_bin[indices2], flux_broad[indices2], linewidth = 0.8, label='Flux broad',color='red')
# plt.plot(wavelength_cube[indices3], th_flux_scaled[indices3], linewidth = 0.5, alpha=0.5, label='Flux no broad')
# plt.xlim(line_lower, line_upper)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Flux [erg/s/cm²/$\AA$]')
# plt.legend(loc='best')
# plt.title('Theoretical vs. Line broadened (Small range)')
# plt.savefig(fig_path+'TheoLineBroadSpec_SmallRange_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='svg')
# plt.show()

# Estimate V-band magnitude
# Multiply spectrum by transmission curve for V band = Filter out all the flux that falls outside the V band
tc = tcs.getTransCurve("Johnson V")
trans = tc(wave_bin)

# Integrate product of spectrum and transmission curve over all wavelengths = Give total flux within the V band
flux_in_v_band = np.trapz(flux_broad * trans, wave_bin)
flux_in_v_band_SED = np.trapz(sed_bin * trans, wave_bin)

# Convert flux in V band to magnitude using zero-point flux density for V band (3.64e-9 erg/s/cm²/A * width A)
v_mag_test = -2.5 * np.log10(flux_in_v_band / flux_0_johnson_v)
v_mag_test_SED = -2.5 * np.log10(flux_in_v_band_SED / flux_0_johnson_v)

print("Estimated V-band magnitude (Test after broadening): ", v_mag_test, v_mag_test_SED)

# Call rebinning for necessary parameters
mirror_reflectivity = rebin_params(wave_color, mirror_ref_nobin, wave_bin, "quadratic")
QE = rebin_params(wave_color, QE_nobin, wave_bin, "quadratic")

# Total transmission (dependent on interpolated mirror reflectivity)
transmission = get_transmission(wave_bin, mirror_reflectivity)

# elemss = []
# for elem in mirror_ref_nobin:
#     elemss.append(elem**nMirrors
#                 * disperser_eff
#                 * cross_disperser_eff
#                 * (1-obscuration)
#                 * FGS_reflectivity
#                 * (1 - pointing_jitter_loss)
#                )

# plt.figure()
# plt.subplot(2, 2, 1)
# plt.plot(wave_bin, mirror_reflectivity, label='interp.')
# plt.scatter(wave_color, mirror_ref_nobin, color='red', s=10, label='colors')
# plt.legend(loc='best')
# plt.title('Interpolated mirror reflectivity')
# # plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Mirror reflectivity')
#
# plt.subplot(2, 2, 2)
# plt.plot(wave_bin, QE, label='interp')
# plt.scatter(wave_color, QE_nobin, color='red', s=10, label='QE')
# plt.legend(loc='best')
# plt.title('Interpolated QE')
# # plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('QE')

# HSA:31/7: not needed cos we have an input flux
# plt.subplot(2, 2, 3)
# plt.plot(wave_bin, IntMag)
# plt.title('Interpolated IntMag')
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('IntMag')

# plt.subplot(2, 2, 4)
# plt.plot(wave_bin, transmission)
# # plt.scatter(wave_color, elemss, color='red', s=10, label='test')
# plt.title('Interpolated transmission')
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Total transmission')
#
# plt.tight_layout()
# plt.show()

# Additionally: straylight baffle from unit photons/s/bin to photons/s/pixel
# straylight_baffle = np.ones(len(wave_bin))       # photons/bin/s (reflected by baffle)  #HSA:31/8
# for wave in range(len(wave_bin)):
#    straylight_baffle[wave] = straylight_baffle[wave] * (resolution/wave_bin[wave]) / 3.08 # Unit: photons/s/pixel

# Blaze function of the spectrograph
blaze_peak = Calc_Blaze_Peak(wave_bin, lambda_cen, lambda_min, lambda_max, order_numbers)

# plt.figure()
# for i in range(len(order_numbers)):
#     wavelength_min = lambda_min[i] * 1e4
#     wavelength_max = lambda_max[i] * 1e4
#
#     plt.axvline(wavelength_min, linestyle="--", linewidth=0.5, color="gray")
#     plt.axvline(wavelength_max, linestyle="--", linewidth=0.5, color="gray")
#
# plt.plot(wave_bin, blaze_peak, linewidth=0.8)
# plt.xlim(lambda_lower, lambda_upper)
# plt.title('Blaze peaks')
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Transmission efficiency')
# plt.savefig('Blaze_peak.svg', facecolor='w', transparent=False, format='svg')
# plt.show()

# Call functions to extract observed simulated spectrum
"""Conversed flux with response (Unit: photons/bin)"""
flux_eff = Calc_ConversionFlux(wave_bin, flux_broad, M1_area, QE, transmission, blaze_peak)
flux, straylight = Calc_Flux_pix(wave_bin, flux_eff, ExpTime, resolution, mean_wavel, transmission, QE)
flux_bin, straylight_bin = Calc_Flux_bin(wave_bin, flux_eff, ExpTime, resolution, mean_wavel, transmission, QE)

sed_eff = Calc_ConversionFlux(wave_bin, sed_bin, M1_area, QE, transmission, blaze_peak)
sed, _ = Calc_Flux_pix(wave_bin, sed_eff, ExpTime, resolution, mean_wavel, transmission, QE)
sed_bin, _ = Calc_Flux_bin(wave_bin, sed_eff, ExpTime, resolution, mean_wavel, transmission, QE)

"""Noise and SNR"""
noise, snr = Calc_Noise(wave_bin, flux, straylight, ExpTime)
noise_bin, snr_bin = Calc_Noise(wave_bin, flux_bin, straylight_bin, ExpTime)
print('Average SNR: ',np.median(snr))
print('Average SNR: ',np.median(snr_bin))

"""Observed simulated spectrum (Units: erg/s/cm²/A) or (photons/bin)"""
ObsSimFlux = Calc_ObsSimFlux(wave_bin, flux, noise, blaze_peak, M1_area, ExpTime, resolution)
ObsSimFlux_bin = Calc_ObsSimFlux(wave_bin, flux_bin, noise_bin, blaze_peak, M1_area, ExpTime, resolution)

# Test new straylight
plt.figure()
# plt.plot(wave_bin, straylight, label='Straylight')
# plt.plot(wave_bin, flux,  label='Flux')
plt.plot(wave_bin, straylight_bin, label='Straylight_bin')
plt.plot(wave_bin, flux_bin,  label='Flux_bin')
plt.xlim(lambda_lower, lambda_upper)
plt.xlabel('Wavelength [$\AA$]')
plt.ylabel('Flux/Straylight [photons/bin]')
plt.legend(loc='best')
# plt.savefig('title.svg', facecolor='w', transparent=False, format='svg')
plt.savefig('title.pdf', facecolor='w', transparent=False, format='pdf')
plt.show()

# Plots
x_ax = np.arange(0, len(wave_bin), 1)

#print(np.mean(np.diff(wave_bin)))

"""Noise"""
plt.figure()
#plt.scatter(x_ax, noise, marker='+', linewidth=0.5, label='Noise')
plt.scatter(x_ax, noise_bin, marker='+', linewidth=0.5, label='Noise_bin', alpha = 0.5)
plt.xlim(x_ax[0], x_ax[-1])
plt.xlabel('Counter')
plt.ylabel('Noise [photons/bin]')
plt.legend(loc='best')
plt.title('Noise')
#plt.savefig(fig_path+'Noise_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='svg')
plt.savefig(fig_path+'Noise_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini),
            facecolor='w', transparent=False, format='pdf')
plt.show()

"""Observed simulated spectrum"""
plt.figure()
#plt.plot(x_ax, ObsSimFlux, linewidth=0.7, label='Flux values')
plt.plot(x_ax, ObsSimFlux_bin, linewidth=0.7, label='Flux')
#plt.scatter(x_ax, sed_bin/blaze_peak, color='orange', label='Abs. Flux', s = 5)
plt.xlim(x_ax[0], x_ax[-1])
plt.xlabel('Counter')
plt.ylabel('Flux [photons/bin]')
plt.legend(loc='best')
plt.title('Observed simulated spectrum')
# plt.savefig(fig_path+'ObservedSimulatedSpectrum_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='svg')
plt.savefig(fig_path+'ObservedSimulatedSpectrum_T{}g{}D{}R{}vsini{}2.pdf'.format(Teff, g, dist_pc, rad_sol, vsini),
            facecolor='w', transparent=False, format='pdf')
plt.show()

"""Signal-to-noise ratio"""
plt.figure()
# plt.plot(wave_bin, snr, linewidth=0.7, label='SNR')
# plt.axhline(y=np.median(snr), color='r', linestyle='--', label='Median SNR '+str(np.round(np.median(snr))))
plt.plot(wave_bin, snr_bin, linewidth=0.7, label='SNR_bin')
plt.axhline(y=np.median(snr_bin), color='r', linestyle='--', label='Median SNR_bin '+str(np.round(np.median(snr_bin))))
plt.xlim(lambda_lower, lambda_upper)
plt.xlabel('Wavelength [$\AA$]')
plt.ylabel('S/N')
plt.legend(loc='best')
plt.title('Signal-to-noise ratio')
# plt.savefig(fig_path+'SNR_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='svg')
plt.savefig(fig_path+'SNR_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini),
            facecolor='w', transparent=False, format='pdf')
plt.show()

"""Normalised spectrum"""
plt.figure()
plt.plot(wave_bin, line_invert(ObsSimFlux_bin/(sed_bin/blaze_peak)), linewidth=0.7, label='Normalised Flux')
plt.xlim(lambda_lower, lambda_upper)
plt.xlabel('Wavelength [$\AA$]')
plt.ylabel('Normalised flux')
#plt.legend(loc='best')
plt.title('Reversed Normalised Flux')
#plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini),
#             facecolor='w', transparent=False, format='svg')
plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini),
            facecolor='w', transparent=False, format='pdf')
plt.show()

# Pulsation time serie
wvl, pulsations = get_pulsations('tutorial_test')        # extract pulsations info
rev_pul = line_invert(pulsations)                        # reverse pulsations
mean_profile = mean_pul(rev_pul)                         # compute mean profile
normed_pul = norm_pul(wvl, rev_pul, mean_profile)        # normalise pulsations

convo_test = signal.fftconvolve(line_invert(ObsSimFlux_bin/(sed_bin/blaze_peak)), normed_pul[:, 56], mode = 'same')

"""Convolution test with a pulsation profile"""
plt.figure()
plt.plot(wave_bin, convo_test, linewidth=0.7, label='Convolved spectrum')
plt.xlim(lambda_lower, lambda_upper)
plt.xlabel('Wavelength [$\AA$]')
plt.ylabel('Normalised flux')
#plt.legend(loc='best')
plt.title('Convo test')
#plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.svg'.format(Teff, g, dist_pc, rad_sol, vsini), facecolor='w', transparent=False, format='svg')
plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini), facecolor='w', transparent=False, format='pdf')
plt.show()

"""Additional test"""
plt.figure()

for i in range(np.shape(normed_pul)[1]):
    convo = signal.fftconvolve(line_invert(ObsSimFlux_bin / (sed_bin / blaze_peak)), normed_pul[:, i], mode='same')
    plt.plot(wave_bin, convo, linewidth = 0.7, alpha = 0.5)
plt.xlim(lambda_lower, lambda_upper)
plt.xlabel('Wavelength [$\AA$]')
plt.ylabel('Normalised flux')
plt.show()

# print(inspect.getsource(pyasl.rotBroad))
# print(inspect.getsource(pyasl.instrBroadGaussFast))