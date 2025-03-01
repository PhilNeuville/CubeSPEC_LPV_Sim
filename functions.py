import os
import shutil

from parameters import *
from PyAstronomy import pyasl
from scipy import signal
from scipy.integrate import simps # integration with Simpson formula
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt

#########################
""" General functions """
#########################


def Get_Tlusty_spec(common_path, final_folder, Teff, g, v):
    """
    Import Tlusty spectrum from file given certain Teff, g, and v

    Args:
    common_path: define the common part of the path to the folder
    final_folder: Addition to .vis depending on which Tlusty spectrum is used (.7 or .17)

    Out:
    Tlusty_spec: 2D array (wvl, flux) of Tlusty spectrum
    """
    """
    Choose final folder for importing Tlusty data (Full = .7-files (all values); SED = .17-files)
    final_folder = 'Grid Full'  # full spectrum
    final_folder = 'Grid SED'  # SED
    """
    if final_folder == 'Grid Full':     # Full spectrum
        filename_add = 7
    elif final_folder == 'Grid SED':    # SED / Absolute flux
        filename_add = 17

    # Construct the full path to the file
    filename = "BG{}g{}v{}.vis.".format(Teff, g, v) + str(filename_add)
    file_path = common_path + '/' + final_folder + '/' + filename

    # Load data from file into a Numpy array
    data = np.loadtxt(file_path)

    # Split the data into 2 columns
    wvl = data[:, 0]                # wavelength [A]
    flx = data[:, 1] * 4 * np.pi    # HSA: converting Eddington flux (H0) to astrophyscial flux

    # Create a 2D array from the two columns
    Tlusty_spec = np.column_stack((wvl, flx))

    return Tlusty_spec


def Get_orders(common_path, order_filename):
    """
    Import the wavelengths information for each spectral order

    Args:
    common_path: define the common part of the path to the folder
    order_filename: name of the .txt file containing info about spectral orders

    Out:
    order_numbers: spectral order numbers
    lambda_cen: central wavelengths
    lambda_min: minimum wavelength
    lambda_max: maximum wavelengths
    """
    # Construct the path of the file
    path_file = common_path + '/' + order_filename

    # Load data from file into a Numpy array
    data = np.loadtxt(path_file, skiprows=1)

    # Split the data into 4 columns
    order_numbers = data[:, 0]      # spectral order numbers
    lambda_cen = data[:, 1]         # central wavelengths [pm = 1e-6m]
    lambda_min = data[:, 2]         # minimum wavelengths [pm = 1e-6m]
    lambda_max = data[:, 3]         # maximum wavelengths [pm = 1e-6m]

    return order_numbers, lambda_cen, lambda_min, lambda_max


def Cut_specrange(wave, flux, wave_min, wave_max):
    """
    Extract the spectral range of CubeSPEC

    Args:
    wave: input wavelength array
    flux: input flux
    wave_min: lower limit of CubeSPEC spectral range
    wave_max: upper limit of CubeSPEC spectral range

    Out:
    wave: CubeSPEC spectral range
    flux: flux corresponding to CubeSPEC spectral range
    """
    # cut the spectrum to CubeSpec wavelength
    spec_range = ((wave >= wave_min) & (wave <= wave_max))
    wave = wave[spec_range]
    flux = flux[spec_range]

    return wave, flux


def Scale_vmag(wave, flux, Vmag):
    """
    Convert flux to certain Vmag

    Args:
    wave [A]: Wavelength array
    flux [erg/s/cm²/A]: Stellar flux array to scale
    Vmag: Wanted V magnitude

    Out:
    wave [A]: Wavelength array
    flux_scale [erg/s/cm²/A]: Theoretical stellar flux corresponding to given Vmag
    """
    # Determine flux density in V band that corresponds to given V magnitude
    v_mag_4_flux_density = flux_0_johnson_v * 10 ** (-0.4 * Vmag)  # in erg/s/cm²

    # Get transmission curve object
    tcs = pyasl.TransmissionCurves()
    # Add passbands
    tcs.addSpitzerIRACPassbands()

    # Determine flux density in spectrum within V band
    tc = tcs.getTransCurve("Johnson V")         # selects Johnson V transmission curve object
    trans = tc(wave)                            # selects transmission curve values [0;1] in wave range
    flux_in_v_band = simps(flux * trans, wave)  # why integrate?: to get the total flux in the V band (from erg/s/cm²/A to erg/s/cm²)

    # Scale spectrum to correspond to V magnitude (of 4)
    scale_factor = v_mag_4_flux_density / flux_in_v_band  # to retrieve the B-Cep mag in V band

    #     flux_scale = []
    #     for w in range(len(wave)):
    #         flux_scale.append(flux[w] * scale_factor)
    #print(scale_factor)
    flux_scale = flux * scale_factor

    return scale_factor, flux_scale # update 13/06
    #return np.array(wave), np.array(flux_scale)


# Rebin the spectrum for broadening (from Julia Bodensteiner)
def Rebin(wave, flux, stepwidth, err=False, verbose=False):
    """
    Function to interpolate from old wavelength (wave) to new wavelength array
    (From Julia Bodensteiner) note: stepwidth was hardcoded as = 1.25
    This is required to use the broadening mechanisms (requiring evenly spaced wavelength array)
    """
    # print('stepwidth = ', stepwidth)
    # wl_rebin with new stepwidth
    if verbose is True:
        print('Rebinning to new stepwidth of %f A.' % stepwidth)

    # define new array based on given stepwidth
    wl_start = wave[0]  # define start
    wl_end = wave[-1]  # and end point of new wave array

    wl_rebin = np.arange(wl_start, wl_end + stepwidth, stepwidth)

    # do the interpolation
    intFunc = interp1d(wave, flux, kind="slinear", fill_value='extrapolate')  # interpolation function
    fl_rebin = np.array(intFunc(wl_rebin))  # rebin the wavelengths (-> evenly spaced)

    if err is not False:
        errFunc = interp1d(wave, err, kind="slinear", fill_value='extrapolate')  # what is err? Use a bool in interp1d?
        err_rebin = np.array(errFunc(wl_rebin))  # is not False: could be a list
        return wl_rebin, fl_rebin, err_rebin

    elif err is False:
        return wl_rebin, fl_rebin


# Line broadening mechanisms: Rotation
def Apply_Rotation(wave, flux, vsini, epsilon=0.6):
    """
    Determine rotational broadening (Convolve original spectrum with this)
    (Rotational broadening, requires evenly spaced wavelength array)

    Args:
    wave [A]: Evenly spaces wavelength array
    flux [erg/s/cm²/A]: Rebinned theoretical stellar flux
    vsini [km/s]: Rotational velocity of star
    epsilon: limb darkening coefficient

    Out:
    flux_rot [erg/s/cm²/A]: Flux corrected for rotational broadening
    """
    flux_rot = pyasl.rotBroad(wave, flux, vsini=vsini, epsilon=epsilon)

    return flux_rot


# Line broadening mechanisms: Instrumental
def Apply_Instrument(wave, flux, resolution):
    """
    Determine instrumental broadening (Convolve already rotational broadened spectrum with this)

    Args:
    wave [A]: Evenly spaces wavelength array
    flux [erg/s/cm²/A]: Theoretical stellar flux (possibly already corrected for rotational broadening)
    resolution: instrument spectral resolution (CubeSPEC ~ 55.000)

    Out:
    flux_instr [erg/s/cm²/A]: Theoretical stellar flux including rotational and instrumental broadening
    """
    flux_instr = pyasl.instrBroadGaussFast(wave, flux, resolution)
    #flux_instr = pyasl.instrBroadGauss(wave, flux, resolution)

    return flux_instr


# Flux-Luminosity-Distance equation scaling
def Calc_TlustyStellarSpec(Tlusty_spec, rad_sol, dist_pc):
    """
    Create stellar spectrum from theoretical Tlusty spectrum

    Args:
    Tlusty_spec: 2D array containing wavelength range and flux [erg/s/cm²/A]
    rad_sol: Radius of star in solar radii
    dist_pc: Distance of star in parsec

    Out:
    wavelength [A]: Wavelength array from Tlusty
    flux_atdist [erg/s/cm²/A]: Theoretical stellar flux
    """
    rad = rad_sol * Rsun_pc  # radius expressed in pc

    # Conversion factor for Tlusty flux to stellar flux at distance ("dilution factor")
    factor = (rad / dist_pc) ** (2)

    # Loop over all flux elements + extract wavelength array
    flux_atdist = []
    wavelength = []
    for i in range(len(Tlusty_spec)):
        flx_tmp = Tlusty_spec[i][1] * factor
        flux_atdist.append(flx_tmp)
        wavelength.append(Tlusty_spec[i][0])

    # PN: Doing the same but without the loop
    # wavelength = Tlusty_spec[:, 0]
    # flux_atdist = Tlusty_spec[:, 1] * factor

    return np.array(wavelength), np.array(flux_atdist)


# New rebinning function for colour dependent parameters in CubeSpec optical characteristics
def Rebin_params(wave, param, wave_bin, kind):
    """
    Interpolate necessary for mirror reflectivity and QE (response detector) -> cos' defined only at several wavelength (PN)

    Args:
    wave [A]: Wavelength array
    param: Parameter array of same length as wave -> The values that need interpolation
    wave_bin [A]: Rebinned wavelength array (from the 'rebin' function defined earlier)
    kind: Type of interpolation (e.g. linear, quadratic, cubic, nearest)

    Out:
    param_interp: Interpolated parameter to new wavelength array (wavebin)
    """
    intFunc = interp1d(wave, param, kind=kind, fill_value='extrapolate')
    param_interp = np.array(intFunc(wave_bin))

    return param_interp

def Get_transmission(wave_bin, mirror_reflectivity):
    """
    Obtain the total transmission (dependent on interpolated mirror reflectivity and QE

    Args:
    mirror_reflectivity: interpolated mirror reflectivity from discrete values

    Out:
    transmission: total transmission
    """
    transmission = []
    for wave in range(len(wave_bin)):
        transmission.append(mirror_reflectivity[wave] ** nMirrors
                            * disperser_eff
                            * cross_disperser_eff
                            * (1 - obscuration)
                            * FGS_reflectivity
                            * (1 - pointing_jitter_loss)
                            )

    #print('Ok it works')
    return transmission

def Calc_Blaze_Peak(wavelength, lambda_cen, lambda_min, lambda_max, order_numbers):
    """
    Calculate the blaze function of the spectrograph

    Args:
    wavelength [A]: Wavelength array
    lambda_cen [pm]: Central wavelength array of orders in spectrograph
    lambda_min [pm]: Lower wavelength boundary array of orders
    lambda_max [pm]: Upper wavelength boundary array of orders
    #ratio_edge_centre_tm: Ratio transmission at edges over transmission in centre

    Out:
    blaze_function: Blaze peaks per order of the spectrograph
    """
    blaze_function = []
    for wave in range(len(wavelength)):
        # Determine order number for given wavelength
        for i in range(len(lambda_min)):
            if lambda_min[i] * 1e4 <= wavelength[wave] and wavelength[wave] <= lambda_max[i] * 1e4:
                m = order_numbers[i]
                index = i

        # Wavelength range of order number
        #wave_min = lambda_min[index] * 1e4  # Unit: Angstrom
        #wave_max = lambda_max[index] * 1e4
        wave_cen = lambda_cen[index] * 1e4

        nu = m * np.pi * (wave_cen - wavelength[wave]) / wavelength[wave]
        blaze_function.append(np.sin(nu) ** 2 / nu ** 2)

    return np.array(blaze_function)


def Calc_ConversionFlux(wave_bin, flux_broad, area, QE, transmission, blaze_peak):
    """
    Convert the theoretical stellar flux from unit [erg/s/cm²/A] to [photons/s/cm²/A]
    Then add the response function

    Args:
    wave_bin [A]: Wavelength array
    flux_broad [erg/s/cm²/A]: Theoretical broadened flux (rotational and instrumental)
    area [cm²]: Detector area
    QE: Quantum efficiency array as function of wavelength
    transmission: Transmission array as function of wavelength
    blaze_peak: Blaze function array containing blaze peaks

    Out:
    flux_eff [photons/s/A]: Effective stellar flux
    """
    FluxCal = []
    for wave in range(len(wave_bin)):
        photon_e = (h * 1e7 * cc) / (wave_bin[wave] * 1e-10)    # Photon energy: h * nu = h * c / lambda [erg/photon] (note: 1J = 1e7 ergs)
        flx_tmp = flux_broad[wave] / photon_e  # Unit: photons/s/cm²/A
        FluxCal.append(flx_tmp)

    # Number of detected photons (includes detector response!): Effective stellar flux
    flux_eff = []
    for wave in range(len(wave_bin)):
        flux_eff.append(FluxCal[wave] * area * QE[wave] * transmission[wave] * blaze_peak[wave])  # Unit: photons/s/A

    return np.array(flux_eff)


def Calc_Flux_pix(wave_bin, flux_eff, ExpTime, resolution, transmission, QE):
    """
    Calculate the amount of photons per pixel that are detected by CubeSPEC and the straylight

    Args:
    wave_bin [A]: Wavelength array
    flux_eff [photons/s/A]: Effective stellar flux
#    IntMag: Magnitude array as function of wavelength  HSA:31/7: not needed cos' we have an input flux
#    SpecType: Spectral type of star
    ExpTime [s]: Exposure time for detection
    resolution [lambda/delta lambda]: Spectral resolution of detector

    Out:
    flux [photons/pixel]: Amount of photons detected
    straylight [photons/pixel]: Amount of straylight photons detected
    """
    # Integrated flux per spectral resolution bin (or resolution element)
    # Spectral width of each bin = wavelength / resolution = dlambda = 5000 / 55000 = 0.091A
    flux = []
    straylight = []
    for wave in range(len(wave_bin)):

        # INFO from GERT email 18/04:
        # The factor "wavelength[color]/resolution” is there to go from photons/nm to photons/bin (or per resolution element).
        # If you want this per pixel instead of per bin, you could just divide by 3.08 (20mu slit width, sampled by 6.5mu pixels).
        """"!!!flux.append(flux_bin * 10**(-0.4 * IntMag[wave]) * ExpTime * nPixels) # Unit: photons/pixel!!!"""

        flux_bin = flux_eff[wave] * wave_bin[wave] / resolution     # Unit: photons/bin/s -> the spectral 'sampling' is wider in the red (definition of a spectrograph)
        flux_pix = flux_bin / nPix_per_resol                        # Unit: photons/pixel/s

        flux.append(flux_pix * ExpTime)  # Unit: photons/pixel #HSA:31/7; already integrated over spatial direction

        # straylight_tmp = straylight_M1 * transmission[wave] * QE[wave] * ExpTime * wave_bin[wave] / resolution / nPix_per_resol # Unit: photons/pixel
        # straylight_tmp2 = straylight_tmp + straylight_baffle * QE[wave] * ExpTime  # !!!! # Unit: photons/pixel
        # straylight.append(straylight_tmp2 + straylight_internal * flux[wave])  # Unit: photons/pixel

        stray_m1 = straylight_M1 * wave_bin[wave] / resolution / nPix_per_resol * transmission[wave] * QE[wave] * ExpTime # Unit: photons/pixel
        stray_baffle = straylight_baffle * QE[wave] * ExpTime                                                             # Unit: photons/pixel
        stray_internal = straylight_internal * flux[wave]                                                                 # Unit: photons/pixel
        straylight.append(stray_m1 + stray_baffle + stray_internal)                                                       # Unit: photons/pixel

    return flux, straylight

def Calc_Flux_bin(wave_bin, flux_eff, ExpTime, resolution, transmission, QE):
    """
    Calculate the amount of photons per bin that are detected by CubeSPEC and the straylight

    Args:
    wave_bin [A]: Wavelength array
    flux_eff [photons/s/A]: Effective stellar flux
    ExpTime [s]: Exposure time for detection
    resolution [lambda/delta lambda]: Spectral resolution of detector
    transmission: Mirror transmission (extrapolated)
    QE: Quantum efficiency (extrapolated)

    Out:
    flux [photons/bin]: Amount of photons detected (per bin)
    straylight [photons/bin]: Amount of straylight photons detected (per bin)
    """
    # Integrated flux per spectral resolution bin (or resolution element)
    # Spectral width of each bin = wavelength / resolution = dlambda = 5000 / 55000 = 0.091A
    flux = []
    straylight = []
    for wave in range(len(wave_bin)):

        # INFO from GERT email 18/04:
        # The factor "wavelength[color]/resolution” is there to go from photons/nm to photons/bin (or per resolution element).
        # If you want this per pixel instead of per bin, you could just divide by 3.08 (20mu slit width, sampled by 6.5mu pixels).
        """"!!!flux.append(flux_bin * 10**(-0.4 * IntMag[wave]) * ExpTime * nPixels) # Unit: photons/pixel!!!"""

        flux_bin = flux_eff[wave] * wave_bin[wave] / resolution     # Unit: photons/bin/s -> the spectral 'sampling' is wider in the red (definition of a spectrograph)

        flux.append(flux_bin * ExpTime)  # Unit: photons/bin #HSA:31/7; already integrated over spatial direction

        stray_m1 = straylight_M1 * wave_bin[wave] / resolution * transmission[wave] * QE[wave] * ExpTime       # Unit: photons/bin
        stray_baffle = straylight_baffle * nPix_per_resol * QE[wave] * ExpTime                                 # Unit: photons/bin
        stray_internal = straylight_internal * flux[wave]                                                      # Unit: photons/bin
        straylight.append(stray_m1 + stray_baffle + stray_internal)                                            # Unit: photons/bin

    return flux, straylight

def Calc_Stray_only(wave_bin, flux, ExpTime, resolution, transmission, QE):
    """
    Calculate the straylight only [photons/bin]

    Args:
    wave_bin [A]: Wavelength array
    flux [photons/bin]:Amount of photons detected
    Exptime [s]: Exposure time for detection
    resolution: Spectral resolution of detector
    transmission: Mirror transmission (extrapolated)
    QE: Quantum efficiency (extrapolated)

    Out:
    straylight [photons/bin]: Amount of straylight photons detected (per bin)
    """
    straylight = []
    for wave in range(len(wave_bin)):
        stray_m1 = straylight_M1 * wave_bin[wave] / resolution * transmission[wave] * QE[wave] * ExpTime  # Unit: photons/bin
        stray_baffle = straylight_baffle * nPix_per_resol * QE[wave] * ExpTime  # Unit: photons/bin
        stray_internal = straylight_internal * flux[wave]  # Unit: photons/bin
        straylight.append(stray_m1 + stray_baffle + stray_internal)

    return straylight


def Calc_Noise(wave_bin, flux, straylight, ExpTime):
    """
    Calculate the noise level in photons/bin given the flux and straylight

    Args:
    wave_bin [A]: Wavelength array
    flux [photons/bin]: Amount of photons detected
    straylight [photons/bin]: Amount of straylight photons detected

    Out:
    noise: Random Gaussian noise based on calculated noise level
    """
    noise = []
    snr = []
    for wave in range(len(wave_bin)):
        noise_level_sqr = flux[wave] + straylight[wave] + nPixels * (read_noise ** 2 + dark_flux * ExpTime) + ((flux[wave] * gain_var) ** 2)
        noise.append(np.random.normal(0, np.sqrt(noise_level_sqr)))
        snr.append(flux[wave] / np.sqrt(noise_level_sqr))

    return np.array(noise), np.array(snr)


# Observed simulated spectrum
def Calc_ObsSimFlux(wave_bin, flux, noise, blaze):
    """
    Function to calculate final observed simulated spectrum in unit [erg/s/cm²/A] [photons/pix]?

    Args:
    wave_bin [A]: Wavelength array
    flux [photons/bin or pixel]: Amount of photons detected
    noise [photons/bin or pixel]: Random Gaussian noise array
    blaze: Blaze function of spectrograph
    area [cm²]: Detector area
    ExpTime [s]: Exposure time of detection
    resolution [lambda/delta lambda]: Spectral resolution of CubeSpec

    Out:
    ObsSimFlux [photons/bin or pixel]: Observed simulated stellar flux as observed by CubeSpec
    """
    ObsSimFlux = []
    for wave in range(len(wave_bin)):
        # flx_tmp = flux[wave]
        ObsSimFlux.append(flux[wave] + noise[wave])
        # ObsSimFlux.append(flux[wave])

    # Correct back for blaze function
    ObsSimFlux = ObsSimFlux / blaze

    return np.array(ObsSimFlux)


def Norm_spec(wave_bin, flux, abs_flux):
    """
    Normalise the spectrum

    Args:
    wave_bin [A]: Wavelength array
    flux [photons/bin or pixel]: Non-normalised flux array
    abs_flux [photons/bin or pixel]: Absolute flux array

    Out:
    normed_flux: Normalised flux
    """
    normed_flux = []
    for wave in range(len(wave_bin)):
        normed_flux.append(flux[wave] / abs_flux[wave])

    return np.array(normed_flux)

#####################################
"""FUNCTIONS FOR PULSATIONS MATTER"""
#####################################

def Get_pulsations(pulsation_dir):
    """
    Extract the pulsational profiles of the time serie delivered by FAMIAS

    Args:
    pulsation_dir: Directory's name of the pulsational time serie of interest 'fx_lm_yz'

    Out:
    wavl: xaxis of the pulsational profiles #TODO: check about to keep RV or turn into wvl!!
    pulsations: ND array (yaxis1, yaxis2, ..., yaxisN) of the whole pulsation time-series data (N files)
    """
    # path to the pulsational time-series (generated with FAMIAS)
    pul_dir_path = pulsation_path + '/' + pulsation_dir

    # # Extract the data of the 'times' file (important for FAMIAS)
    # time_filename = "times.dat"
    # time_file_path = pul_dir_path + '/' + time_filename
    # # Load data from file into a Numpy array
    # time_data = np.genfromtxt(time_file_path)
    # # Split the data into columns
    # # is 3rd column the weight of each spectrum? #TODO (PN: I think so)
    # files_id, time_span, weights = time_data[:, 0], time_data[:, 1], time_data[:, 2]

    # Count number of files in the directory (to loop automatically over the time-series)
    count = 0
    for path in os.listdir(pul_dir_path):
        if os.path.isfile(os.path.join(pul_dir_path, path)):
            count += 1

    #print('Number of files in folder "{}": {} files'.format(pulsation_dir, count))

    # Loop over pulsational time-series
    for i in range(count - 1):
        # Construct the full path to the file
        filename = "{}.dat".format(i)
        file_path = pul_dir_path + '/' + filename
        #print(filename)
        # Load data from file into a Numpy array
        data = np.loadtxt(file_path)

        if i == 0:
            # Split the data into 2 columns
            rv = data[:, 0]        # x_axis
            pulsations = data[:, 1]      # y_axis: normalised flux
            # Create an (2D) array from the columns
            #pulsations = np.column_stack((xaxis, yaxis))

        else:
            # Data of the i-th pulsation file
            yaxis = data[:, 1]      # y_axis: normalised flux
            # Create an (ND) array by adding a column
            pulsations = np.column_stack((pulsations, yaxis))

    # return np.array(wvl), pulsations
    return rv, pulsations


def Get_null_profile(null_dir):
    """
    Extract the null pulsation profiles implemented in FAMIAS (constant gaussian profile)

    Args:
    null_dir: Directory's name of the null profile pulsation

    Out:
    wavl: xaxis of the pulsational profiles TODO: check about to keep RV or turn into wvl
    pulsations: Np array of the null pulsation profile
    """
    # path to the pulsational time-series (from FAMIAS)
    null_dir_path = pulsation_path + '/' + null_dir

    # Construct the full path to the file
    filename = "0.dat" # Need only one file (all the same in the null time serie)
    file_path = null_dir_path + '/' + filename

    # Load data from file into a Numpy array
    data = np.loadtxt(file_path)

    # Split the data into 2 columns
    wvl = data[:, 0]  # x_axis
    null_pulsations = data[:, 1]  # y_axis: normalised flux

    return np.array(wvl), np.array(null_pulsations)


def Line_invert(spectral_array):
    """
    Invert a normalised line/spectrum (ones to zeroes/zeroes to ones)

    Args:
    spectral_array: Array containing a spectrum, a line or a pulsational profile

    Out:
    reversed_array: Reversed line/spectrum/profile
    """
    # Subtract the spectral array from a similar made of ones
    reversed_array = np.ones_like(spectral_array) - spectral_array

    return reversed_array


def Mean_pul(pulsations):
    """
    Compute the mean profile of the time serie
    ((and use it for deconvolving the time serie (to remove pressure broadening?)))

    Args:
    pulsations: ND array (yaxis1, yaxis2, ..., yaxisN) of the whole pulsation time series data (N files)

    Out:
    mean_profile: Mean profile of the pulsation profiles time serie
    """
    mean_profile = []
    # Compute the average profile of the pulsational time serie
    for i in range(np.shape(pulsations)[0]):
        mean_pt = np.mean(pulsations[i, :])
        mean_profile.append(mean_pt)

    # Deconvolve each time series step with the mean profile
    # for j in range(np.shape(pulsations)[1]):
    #     pulsations[:, j] = signal.deconvolve(pulsations[:, j], mean_profile)

    return np.array(mean_profile)


def Norm_pul(wvl, pul, mean_profile):
    """
    Normalise the pulsations (making a so-called "kernel") #TODO: use the mean profile or each step's one?

    Args:
    wvl: xaxis of the pulsation profiles
    pulsations: ND array (yaxis1, yaxis2, ..., yaxisN) of the whole pulsation time series data (N files)
    mean_profile: Mean pulsation profile of the time serie

    Out:
    wvl: xaxis of the pulsational time series
    norm_pulsations: Normalised pulsations (still ND array)
    """
    # New variable to work with
    normed_pulsations = np.zeros_like(pul)

    # Integrate the mean pulsation profile
    mean_int = np.trapz(mean_profile, wvl)

    for i in range(np.shape(pul)[1]):
        # Integrate each pulsational profile and divide by it
        #pul_int = np.trapz(normed_pulsations[:, i], wvl)
        #normed_pulsations[:, i] = normed_pulsations[:, i] / pul_int
        # Normalise by the mean profile
        normed_pulsations[:, i] = pul[:, i] / mean_int

    return normed_pulsations

def Convolve_spec_puls(normalised_spectrum, puls, mode):
    """
    Make the convolution between the (normalised & inverted) Tlusty spectrum and (normalised & inverted) pulsation profiles

    Args:
    normalised_spectrum: Normalised, inverted, and scaled Tlusty spectrum
    puls: Normalised and inverted pulsation profiles (ND array)
    mode: string indicating the size of the output ('full', 'valid' or 'same')

    Out:
    convolved_spectra: Time serie of pulsating spectra
    """
    # New ND array to work with
    convolved_spectra = np.zeros((len(normalised_spectrum), np.shape(puls)[1]))

    # Convolving the spectrum with each pulsation step
    for pul in range(np.shape(puls)[1]):
        print(pul + 1, ' / ', np.shape(puls)[1])
        #convolved_spectra[:, pul] = signal.fftconvolve(normalised_spectrum, puls[:, pul], mode)
        convolved_spectra[:, pul] = Line_invert(signal.fftconvolve(normalised_spectrum, puls[:, pul], mode))
        #convolved_spectra[:, pul] = Line_invert(convolved_spectra[:, pul])

    return convolved_spectra


def Convolve_spec_puls_log(rv, wave_bin, normalised_spectrum, puls, mode):
    """
    Make the convolution between the (normalised & inverted) Tlusty spectrum and (normalised & inverted) pulsation profiles
    Note: THIS IMPLEMENTS THE METHOD OF TIMOTHY INVOLVING ln(wvl)
    Args:
    rv: Radial velocity vector of the pulsation profile [km/s]
    wave_bin: Wavelength array [A]
    normalised_spectrum: Normalised, inverted, and scaled Tlusty spectrum
    puls: Normalised and inverted pulsation profiles (ND array)
    mode: string indicating the size of the output ('full', 'valid' or 'same')

    Out:
    convolved_spectra: Time serie of pulsating spectra
    """
    # Move to the logspace of wvl
    wave_log = np.log(wave_bin)
    # ln(wvl) = int(v/c)
    rv_c = rv / (cc * 1e-3)
    # rebin the log(lambda) with the same stepsize as the v/c array
    wave_log_bin, spec_log_bin = Rebin(wave_log, normalised_spectrum, np.diff(rv_c)[0])
    # Make the convolution
    convolution_tmp = Convolve_spec_puls(spec_log_bin, puls, mode)
    # Back to wvl space
    wave_back = np.exp(wave_log_bin)
    # Rebin the convolved spectra back on the initial wavelength array
    for step in range(np.shape(convolution_tmp)[1]):
        print(step)
        if step == 0:
            _, convolved_spectra = Rebin(wave_back, convolution_tmp[:, step], mean_wavel/resolution/nPix_per_resol)
            # We remove the last element as the Rebin function adds one
            # This is to keep using the wave_bin vector as done from the start
            convolved_spectra = convolved_spectra[:-1]
        else:
            _, convolved_spectrum_step = Rebin(wave_back, convolution_tmp[:, step], mean_wavel/resolution/nPix_per_resol)
            # We remove the last element as the Rebin function adds one
            convolved_spectrum_step = convolved_spectrum_step[:-1]
            convolved_spectra = np.column_stack((convolved_spectra, convolved_spectrum_step))

    return convolved_spectra


def Rebin_puls(rv, normed_puls, line_wvl):
    """
    Take the normalised pulsation kernel and returns it in the wavelength space
    Rebinned with the same stepsize as the one used for rebinning the Tlusty spectrum
    stepsize = mean_wavel / resolution / nPix_per_resol

    Args:
    rv: Xaxis of pulsation profiles [km/s]
    normed_puls: Normalised pulsation profiles [norm. flux]
    line_wvl: Wavelength of the spectral line of interest

    Out:
    wvl_bin: Rebinned wavelength xaxis of pulsational profiles
    pul_bin: Rebinned ND array of the pulsation profiles
    """
    # Wvl array of pulsation profile acc to the line
    wvl = (rv / (cc * 1e-3) + 1) * line_wvl

    for i in range(np.shape(normed_puls)[1]):
        if i == 0:
            # Rebinning each pulsation with the same stepsize
            wvl_bin, pul_bin = Rebin(wvl, normed_puls[:, i], mean_wavel/resolution/nPix_per_resol)

        else:
            _, var = Rebin(wvl, normed_puls[:, i], mean_wavel/resolution/nPix_per_resol)
            # Create an (ND) array by adding a column
            pul_bin = np.column_stack((pul_bin, var))

    return wvl_bin, pul_bin


def Select_lines(full_xaxis, full_yaxis, lower_lim, upper_lim):
    """
    Take the entire spctra time series and split it in line(s)

    Args:
    full_xaxis: Entire xaxis quantity [A]
    full_yaxis: Entire yaxis quantity [Normalised flux]
    lower_lim: List of lower spectral line limits [A]
    upper_lim: List of upper spectral line limits [A]

    Out:
    To complete...
    """
    return

########################
""" Deal with outputs"""
########################


def Give_outputs(output_path, folder_name, wvl, pulsations, wvl_inf, wvl_sup):
    """
    Write the (output) data into '.dat' files in a given folder

    Args:
    output_path: Path to the output folder
    folder_name: Name of the folder to write the output files in
    wvl: Data xaxis quantity - wavelength array [A]
    pulsations: Data yaxis quantity - pulsating spectra [norm. flux]
    wvl_inf: Inferior bound of the line spectral range [A]
    wvl_sup: Superior bound of the line spectral range [A]

    Out:
    /
    """
    # Create the folder for the output files
    out_path_dir = output_path + '/' + folder_name + '/'
    # Check whether the specified path exists or not
    isExist = os.path.exists(out_path_dir)
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(out_path_dir)
        print("The new OUTPUT directory: " + folder_name + ", is created!")

    # # Store the data by stacking the axes in 2 columns
    # DataOut = np.column_stack((xaxis, yaxis))
    # # Write the data into a '.dat' file
    # np.savetxt(out_path_dir + 'test.dat', DataOut)

    # Copy the times.dat in the output directory (accounts for file names and cadence [days])
    shutil.copyfile(pulsation_path + '/' + pulsation_dir + '/' + 'times.dat', out_path_dir + 'times.dat')

    # Extract the spectral range of the line
    spec_range = ((wvl >= wvl_inf) & (wvl <= wvl_sup))
    wvl = wvl[spec_range]

    # Works for data that are stored
    for i in range(np.shape(pulsations)[1]):
        # Store the data by stacking the axes in 2 columns
        DataOut = np.column_stack((wvl, pulsations[:, i][spec_range]))
        # Write the data into a '.dat' file
        np.savetxt(out_path_dir + '{}.dat'.format(i), DataOut, "%.8f")

    return


def Give_times(output_path, folder_name, pulsation_dir):
    """
    Write the times.dat file required by FAMIAS to call the time-series
    I'd like to make it flexible to indtoruce interruptions in the data.

    Args:
    output_path: Path to the output folder
    folder_name: Name of the folder to write the times.dat file in
    pulsation_dir: Folder containing the time series delivered by FAMIAS

    Out: /
    """
    # Create the folder for the output files
    out_path_dir = output_path + '/' + folder_name + '/'

    # Copy the times.dat in the output directory (accounts for file names and cadence [days])
    shutil.copyfile(pulsation_path + '/' + pulsation_dir + '/' + 'times.dat', out_path_dir + 'times.dat')

    return


"""MISCELLANEOUS"""


def Huge_all(wave_bin, pulsation_dir, rv, normalised_spectrum, puls, mode, sed_bin, transmission, QE, blaze_peak, line_names, line_infs, line_sups):

    # First create the folder and copy/paste the times.dat file in it
    for line1 in range(len(line_names)):
        # Create the folder for the output files
        out_path_dir = output_path + '/' + pulsation_dir + '/' + line_names[line1] + '/'
        # Check whether the specified path exists or not
        isExist1 = os.path.exists(out_path_dir + 'noise' + '/')
        isExist2 = os.path.exists(out_path_dir + 'noiseless' + '/')
        if not isExist1:
            # Create a new directory because it does not exist
            os.makedirs(out_path_dir + 'noise' + '/')
            # print("The new OUTPUT directory: " + line_names[line1] + ", is created!")
        if not isExist2:
            # Create a new directory because it does not exist
            os.makedirs(out_path_dir + 'noiseless' + '/')
            # print("The new OUTPUT directory: " + line_names[line1] + ", is created!")

        # Copy the times.dat in the output directory (accounts for file names and cadence [days])
        shutil.copyfile(pulsation_path + '/' + pulsation_dir + '/' + 'times.dat', out_path_dir + 'noise' + '/' + 'times.dat')
        shutil.copyfile(pulsation_path + '/' + pulsation_dir + '/' + 'times.dat', out_path_dir + 'noiseless' + '/' + 'times.dat')
        print('Folders for', line_names[line1], ' and times.dat file: DONE')

    # Move to the logspace of wvl
    wave_log = np.log(wave_bin)
    # ln(wvl) = int(v/c)
    rv_c = rv / (cc * 1e-3)
    # Rebin the log(wvl) with the same stepsize as the v/c array
    wave_log_bin, spec_log_bin = Rebin(wave_log, normalised_spectrum, np.diff(rv_c)[0])
    # Back to wvl space
    wave_back = np.exp(wave_log_bin)
    for pul in range(np.shape(puls)[1]):
        # Make the convolution btw TLUSTY and a pulsation step
        convolution_tmp = Line_invert(signal.fftconvolve(spec_log_bin, puls[:, pul], mode))
        # Rebin the spectrum with the usual wvl_bin
        _, convolved_spectrum = Rebin(wave_back, convolution_tmp, mean_wavel/resolution/nPix_per_resol)
        # Remove last element as the Rebin function adds one (this is to keep using wave_bin vector)
        convolved_spectrum = convolved_spectrum[:-1]

        # if noise_status == True:
        # print('Step nb', pul + 1, ' / ', np.shape(puls)[1], ' NOISE START')
        # Move back to physical flux to compute the photon flux required to compute the noise
        phys_flux = convolved_spectrum * sed_bin
        # Compute straylight
        straylight = Calc_Stray_only(wave_bin, phys_flux, ExpTime, resolution, transmission, QE)
        # Use the flux and straylight to compute the noise distribution & SNR
        noise, snr = Calc_Noise(wave_bin, phys_flux, straylight, ExpTime)
        # Combine the actual signal with the noise distribution to obtain the observed signal
        obsSimFlux = Calc_ObsSimFlux(wave_bin, phys_flux, noise, blaze_peak)
        # Normalise the observed flux
        obsSimFlux_normed = Norm_spec(wave_bin, obsSimFlux, sed_bin / blaze_peak)
        # Rounding
        final_flux_noise = np.round(obsSimFlux_normed, 8)
        # print('Step nb', pul + 1, ' / ', np.shape(puls)[1], ' NOISE DONE')
        # plt.plot(wave_bin, final_flux, lw = 0.7)
        # plt.show()

        # else:
        # print('Step nb', pul + 1, ' / ', np.shape(puls)[1], ' NOISELESS START')
        final_flux_noiseless = np.round(convolved_spectrum, 8)
        # print('Step nb', pul + 1, ' / ', np.shape(puls)[1], ' NOISELESS DONE')
        # plt.plot(wave_bin, final_flux, lw = 0.7)
        # plt.show()

        for line in range(len(line_names)):

            # Extract the spectral range of the line
            spec_range = ((wave_bin >= line_infs[line]) & (wave_bin <= line_sups[line]))
            # wave_bin = wave_bin[spec_range]

            # Stack the data in columns
            DataOut_noise = np.column_stack((wave_bin[spec_range], final_flux_noise[spec_range]))
            DataOut_noiseless = np.column_stack((wave_bin[spec_range], final_flux_noiseless[spec_range]))

            # Redefine the output directory path
            out_path_dir = output_path + '/' + pulsation_dir + '/' + line_names[line] + '/'
            # Write the data into a '.dat' file
            np.savetxt(out_path_dir + 'noise' + '/' + '{}.dat'.format(pul), DataOut_noise, "%.8f")
            np.savetxt(out_path_dir + 'noiseless' + '/' + '{}.dat'.format(pul), DataOut_noiseless, "%.8f")

            # plt.plot(wave_bin[spec_range], final_flux[spec_range], lw=0.7)
            # plt.show()

        print('Step nb', pul + 1, ' / ', np.shape(puls)[1], ' DONE')

    return