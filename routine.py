# Imports
from functions import *
from parameters import *
import matplotlib.pyplot as plt

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

"""Tlusty spectrum scaling and normalisation"""

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

# Call broadening mechanisms
flux_rot = Apply_Rotation(wave_bin, flux_bin, vsini) # Rotation
flux_broad = Apply_Instrument(wave_bin, flux_rot, resolution) # Instrument
#TODO: il y a des effets de bords sur le spectre normalisé?! -> ça semble due au instrumental broadening

# Compute the blaze function of the spectrograph
blaze_peak = Calc_Blaze_Peak(wave_bin, lambda_cen, lambda_min, lambda_max, order_numbers)

# Call rebinning & interpolation for necessary parameters (reflectivity, quantum efficiency and transmission)
mirror_reflectivity = Rebin_params(wave_color, mirror_ref_nobin, wave_bin, "quadratic")
QE = Rebin_params(wave_color, QE_nobin, wave_bin, "quadratic")
transmission = Get_transmission(wave_bin, mirror_reflectivity)

# Extract effective flux/SED [photons/A/s]
flux_eff = Calc_ConversionFlux(wave_bin, flux_broad, M1_area, QE, transmission, blaze_peak)
sed_eff = Calc_ConversionFlux(wave_bin, sed_bin, M1_area, QE, transmission, blaze_peak)
#TODO: watch out that the normalisation gets rid of the blaze_peak shape!
# To keep in mind when returning the observed spectrum (signal + noise)/blaze_peak

# Computes the flux and nominal flux (SED) per spectral bin [photons/bin]
flux_bin, _ = Calc_Flux_bin(wave_bin, flux_eff, ExpTime, resolution, mean_wavel, transmission, QE)
sed_bin, _ = Calc_Flux_bin(wave_bin, sed_eff, ExpTime, resolution, mean_wavel, transmission, QE)

# Normalise the spectrum
normed_spec = Norm_spec(wave_bin, flux_bin, sed_bin)
rev_normed_spec = Line_invert(normed_spec)

# Size of the plot figures
plt.rcParams['figure.figsize'] = [12, 5]

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

"""Pulsation time serie delivered by FAMIAS: inversion and normalisation"""

# Import the pulsation profiles
rv, puls = Get_pulsations(pulsation_dir)

# Invert pulsation profiles
rev_puls = Line_invert(puls)

# Compute mean (inverted) pulsation profile
mean_pul_profile = Mean_pul(rev_puls)

# Normalise pulsations (w/ mean or static profile?)
normed_puls = Norm_pul(rv, rev_puls, mean_pul_profile)

### Plots the (normalised) (inverted) pulsation profiles
# plt.figure()
# for i in range(np.shape(puls)[1]):
#     plt.plot(rv, normed_puls[:, i], alpha = 0.2)
# plt.xlabel('RV [km/s]')
# plt.ylabel('Normalised flux')
# #plt.legend(loc='best')
# plt.title('Pulsation profiles')
# plt.show()

""" Convolution between the normalised spectrum and the pulsation profiles"""

# Call the convolution between the spectrum and the pulsation profiles (ND array)
convolved_spectra = Convolve_spec_puls(rev_normed_spec, normed_puls, 'same')

plt.figure()
# plt.plot(wave_bin, normed_spec + 0.2, linewidth=0.7, label = 'Normalised spectrum')
for pul in range(np.shape(convolved_spectra)[1]):
    plt.plot(wave_bin, convolved_spectra[:, pul], linewidth=0.7, label = 'Convolved spectrum')
# plt.axvline(4471.292, lw = 0.5, ls = '--')
# plt.axvline(4861.3, lw = 0.5, ls = '--')
# plt.legend(loc='best')
plt.show()

""" Computes the noise & SNR for each spectrum to obtain the Observed spectrum"""

#TODO: obs = (signal + noise) / blaze_peak

""" Write the time series of pulsing spectra as outputs in .dat files"""

# Call the function to save each pulsting spectrum step as .dat files in a given folder
Give_outputs(output_path, 'routine_test', wave_bin, convolved_spectra)