# Imports
import imageio.v2 as imageio

from datetime import datetime
from functions import *
from parameters import *
import matplotlib.pyplot as plt

# Start the clock
start_time = datetime.now()

# Size of the plot figures
plt.rcParams['figure.figsize'] = [12, 5]

# Create a new "figures" folder for the star of interest
isExist = os.path.exists(fig_path)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(fig_path)
    print("The new FIGURES directory is created!")

# Create a new "output" folder for the star of interest
isExist_out = os.path.exists(output_path)
if not isExist_out:
    # Create a new directory because it does not exist
    os.makedirs(output_path)
    print("The new OUTPUT directory is created!")

"""
Full routine:
Make use of the Tlusty spectrum and FAMIAS pulsation time-serie
Convolve the 'treated' spectrum with pulsations

Inputs:
Tlusty spectrum
Parameters
FAMIAS pulsation time serie

Output: time serie of spectra affected by pulsations and written in separated files
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
scale_factor, th_flux_scaled = Scale_vmag(wavelength_cube, th_flux_cube, Vmag)
th_SED_scaled = th_SED_cube * scale_factor                                      # Update 13/06
# There was a PROBLEM somewhere, the SED was smaller than the flux after the scaling (due to different flux in V band)
# the scaling factors were different and so the scaled SED was smaller than the scaled flux
# the normalisation resulted in a continuum > 1
# _, th_SED_scaled = Scale_vmag(wavelength_cube, th_SED_cube, Vmag)

# Call rebinning function (using the spectrograph )
wave_bin, flux_bin = Rebin(wavelength_cube, th_flux_scaled, mean_wavel/resolution/nPix_per_resol)
_, sed_bin = Rebin(wavelength_cube, th_SED_scaled, mean_wavel/resolution/nPix_per_resol)

"""
# OLD LOCATION: Call broadening mechanisms (to be commented out if done already through FAMIAS)
# flux_rot = Apply_Rotation(wave_bin, flux_bin, vsini) # Rotation
# flux_broad = Apply_Instrument(wave_bin, flux_rot, resolution) # Instrument
#TODO: edge effects on normalised spectrum?! -> due to instrumental broadening (PyAstronomy)
"""

# Compute the blaze function of the spectrograph
blaze_peak = Calc_Blaze_Peak(wave_bin, lambda_cen, lambda_min, lambda_max, order_numbers)

# Call rebinning & interpolation for necessary parameters (reflectivity, quantum efficiency and transmission)
mirror_reflectivity = Rebin_params(wave_color, mirror_ref_nobin, wave_bin, "quadratic")
QE = Rebin_params(wave_color, QE_nobin, wave_bin, "quadratic")
transmission = Get_transmission(wave_bin, mirror_reflectivity)

# Extract effective flux/SED [photons/A/s]
flux_eff = Calc_ConversionFlux(wave_bin, flux_bin, M1_area, QE, transmission, blaze_peak)
sed_eff = Calc_ConversionFlux(wave_bin, sed_bin, M1_area, QE, transmission, blaze_peak)
""" Watch out that the normalisation by the SED gets rid of the blaze_peak shape! """
""" To keep in mind when returning the observed spectrum (signal + noise)/blaze_peak """

# Computes the flux and nominal flux (SED) per spectral bin [photons/bin]
flux_bins, _ = Calc_Flux_bin(wave_bin, flux_eff, ExpTime, resolution, transmission, QE)
sed_bin, _ = Calc_Flux_bin(wave_bin, sed_eff, ExpTime, resolution, transmission, QE)

### Plots the not-(yet)-normalised spectrum
# plt.figure()
# plt.plot(wave_bin, flux_bins, linewidth=0.7, label='Flux')
# plt.plot(wave_bin, sed_bin, linewidth=0.7, label='SED')
# plt.xlim(lambda_lower - 10, lambda_upper + 10)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Normalised flux')
# plt.legend(loc='best')
# plt.title('Normalised Flux')
# # plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini), facecolor='w', transparent=False, format='pdf')
# plt.show()

# Normalise the spectrum
normed_spec = Norm_spec(wave_bin, flux_bins, sed_bin)

# Invert the normalised spectrum
rev_normed_spec = Line_invert(normed_spec)

# NEW LOCATION: Call broadening mechanisms (to be commented out if done already through FAMIAS)
# flux_rot = Apply_Rotation(wave_bin, normed_spec, vsini) # Rotation
# flux_broad = Apply_Instrument(wave_bin, flux_rot, resolution) # Instrument
# rev_normed_spec = Line_invert(flux_rot)

### Plots the normalised spectrum
plt.figure()
plt.plot(wave_bin, normed_spec, linewidth=0.7)
# plt.scatter(wave_bin, normed_spec, s = 5)
plt.axhline(1, ls = '--', c='k')
plt.xlim(lambda_lower, lambda_upper)
plt.xlabel('Wavelength [$\AA$]')
plt.ylabel('Normalised flux')
plt.title('Normalised Flux')
# plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini), facecolor='w', transparent=False, format='pdf')
plt.close()

print('Tlusty scaling and normalisation: DONE \n')

###########################################################################
"""Pulsation time serie delivered by FAMIAS: inversion and normalisation"""
###########################################################################


# Import the pulsation profiles
# rv, puls = Get_pulsations(pulsation_dir)
rv, puls = Get_pulsations(dir)

# Invert pulsation profiles
rev_puls = Line_invert(puls)

# Compute mean (inverted) pulsation profile
mean_pul_profile = Mean_pul(rev_puls)

# Normalise pulsations (w/ mean or static profile?)
normed_puls = Norm_pul(rv, rev_puls, mean_pul_profile)

print('Pulsation time series: DONE \n')

#############################################################################
""" Convolution between the normalised spectrum and the pulsation profiles"""
#############################################################################

# Call the convolution between spectrum and pulsation profiles (ND array)
# convolved_spectra = Convolve_spec_puls(rev_normed_spec, normed_puls, 'same') # Initial one

# Switch RV to wavelength and rebin the pulsation profiles with the same spectral stepsize
# wvl_bin, pul_bin = Rebin_puls(rv, normed_puls, line_wvl)
# convolved_spectra2 = Convolve_spec_puls(rev_normed_spec, pul_bin, 'same')    # By attempting to rebin in the wvl space

# Make the convolution between spectrum and pulstaion profiles IN THE log(wvl) SPACE (TIMOTHY METHOD)
convolved_spectra_log = Convolve_spec_puls_log(rv, wave_bin, rev_normed_spec, normed_puls, 'same')

### Plots the normalised spectrum
# plt.figure()
# plt.plot(wave_bin, normed_spec, linewidth=0.7, label='Normalised Flux', alpha = 0.5)
# # plt.plot(wave_bin, convolved_spectra[:, 0], linewidth = 0.7, label='Conv. in RV space')
# # plt.plot(wave_bin, convolved_spectra2[:, 0], linewidth = 0.7, label='Conv. after rebin in wvl space', c ='r')
# # plt.plot(wave_bin, convolved_spectra_log[:, 0], linewidth = 0.7, label='Conv. Timothy function', ls ='--')
# plt.plot(wave_bin, convolved_spectra_log[:, 0], linewidth = 0.7, c = 'r')
# plt.xlim(lambda_lower, lambda_upper)
# plt.xlabel('Wavelength [$\AA$]')
# plt.ylabel('Normalised flux')
# plt.legend(loc='best')
# plt.title('Normalised Flux')
# for i in range(len(line_name)):
#     plt.axvline(line_cen[i], lw = 0.5, ls = '--', c = 'grey')
#     plt.text(line_cen[i], 1.02, line_name[i], fontsize=12)
# # plt.savefig(fig_path+'normFlux_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, dist_pc, rad_sol, vsini), facecolor='w', transparent=False, format='pdf')
# plt.close()

print('Convolutions |b| Tlusty & pulsations: DONE \n')

#################################################################################
""" Computes the noise & SNR for each spectrum to obtain the Observed spectrum"""
#################################################################################

# New ND array to store the convolved, noisy, normalised spectra to be later written in files
final_time_series = np.zeros_like(convolved_spectra_log)
# The noise is random, so it needs to be re-computed for each step of the time-series
for step in range(np.shape(final_time_series)[1]):

# However, the noise requires to use the photon flux! So we need to move back to the physical flux again
    phys_flux_step = convolved_spectra_log[:, step] * sed_bin

# Compute the straylight for each spectrum of the time-series
    straylight_step = Calc_Stray_only(wave_bin, phys_flux_step, ExpTime, resolution, transmission, QE)

# Use the flux and straylight to compute the noise distribution
    noise_step, snr_step = Calc_Noise(wave_bin, phys_flux_step, straylight_step, ExpTime)

# Combine the actual signal with the noise distribution to obtain the observed signal
    ObsSimFlux = Calc_ObsSimFlux(wave_bin, phys_flux_step, noise_step, blaze_peak)

# Normalise the observed flux
    ObsSimFlux_normed = Norm_spec(wave_bin, ObsSimFlux, sed_bin / blaze_peak)

# Storing the step spectrum in an ND array accounting for the whole time series
    final_time_series[:, step] = np.round(ObsSimFlux_normed, 8)

print('Noise computation: DONE \n')

############################################################
# Loop over the list of line to save (figs, Gifs and data) #
############################################################
for line in range(len(line_name)):

    ########################################################################
    """ Write the time series of pulsing spectra as outputs in .dat files"""
    ########################################################################

    # Call the function to save each pulsating spectrum step as .dat files in a given folder
    Give_outputs(output_path, line_name[line], wave_bin, final_time_series, line_inf[line], line_sup[line])

    # Look at the clock
    # mid_time = datetime.now()
    # # Current running time
    # print('Runtime:', mid_time - start_time)

    #########################################################################
    """Save the plots of the time series and make Gifs animation out of it"""
    #########################################################################

    ## Put all of this in a function?
    # folder_name = 'log_vsini100'
    folder_name = line_name[line]

    # Create the folder for the output files
    fig_path_dir = fig_path + '/' + 'Animations' + '/' + folder_name
    # Check whether the specified path exists or not
    isExist = os.path.exists(fig_path_dir)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(fig_path_dir)
        print("The new FIGURES directory: " + folder_name + ", is created!")

    #spec_range = ((wave_bin >= line_inf[line]) & (wave_bin <= line_sup[line]))
    # Plot the LPV and the pulsational profiles
    filenames_mix = []
    for i in range(np.shape(final_time_series)[1]):
        fig, ax = plt.subplots(2, 1, dpi = 150, figsize = (8, 10))
        ax[0].plot(wave_bin, final_time_series[:, i], linewidth=0.7)
        ax[0].set_xlim(line_inf[line], line_sup[line])
        ax[0].set_ylim(0.95 * np.min(final_time_series[:, 0][((wave_bin >= line_inf[line]) & (wave_bin <= line_sup[line]))]), 1.04)
        ax[0].set_xlabel('Wavelength [$\AA$]')
        ax[0].set_ylabel('Normalised flux')
        ax[0].set_title('LPV ' + line_name[line])
        #
        ax[1].plot(rv, puls[:, i], linewidth=0.7)
        ax[1].set_xlabel('RV [km/s]')
        ax[1].set_ylabel('Normalised flux')
        ax[1].set_ylim(0.98 * np.min(puls[:, 0]), 1.01)
        ax[1].set_title('Pulsation profiles')
        #
        plt.savefig(fig_path_dir + f'/{i}.png', facecolor='w', transparent=False, dpi=150)
        filename_mix = fig_path_dir + f'/{i}.png'
        filenames_mix.append(filename_mix)
        # plt.tight_layout()
        plt.close()

    ## Creating the GIFs from previously saved figures
    with imageio.get_writer(fig_path + '/' + 'Animations' + '/' + '{}.gif'.format(folder_name), fps=20) as writer:
        for filename in filenames_mix:
            image = imageio.imread(filename)
            writer.append_data(image)
    print(line_name[line] + ' line GIF saved!')

    shutil.rmtree(fig_path_dir) # Delete the figures used for the GIF (save space)

# Look at the clock
end_time = datetime.now()
# Full running time
print('Runtime:', end_time - start_time)