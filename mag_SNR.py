# Imports
from datetime import datetime
from functions import *
import matplotlib.pyplot as plt

"""
Compute the SNR for different magnitudes in the V band
"""

# Start the clock
start_time = datetime.now()

# Import Tlusty spectrum from file given certain Teff, g, v
final_folder = 'Grid Full'
Tlusty_spec = Get_Tlusty_spec(common_path, final_folder, Teff, g, v)
# Import the wavelengths for each order
order_numbers, lambda_cen, lambda_min, lambda_max = Get_orders(common_path, order_filename)

# Call theoretical spectrum (scaled with dilution factor)
wavelength, th_flux = Calc_TlustyStellarSpec(Tlusty_spec, rad_sol, dist_pc)

# Select only CubeSPEC wavelength range & corresponding flux
wavelength_cube, th_flux_cube = Cut_specrange(wavelength, th_flux, lambda_lower, lambda_upper)

# V-mag array
step = 1
min = 0
max = 8
mag_list = np.arange(min, max + step, step)
exptime_list = [60, 300, 900]

# Plot SNR vs V-mag
plt.figure()

for ExpTime in exptime_list:
    SNR_ls = []
    for Vmag in mag_list:
        print(Vmag)

        # Scale Tlusty spectrum to correspond to the V magnitude given as input
        _, th_flux_scaled = Scale_vmag(wavelength_cube, th_flux_cube, Vmag)

        # Call rebinnning function (required for broadeing mechanism convolutions)
        wave_bin, flux_bin = Rebin(wavelength_cube, th_flux_scaled, stepwidth = mean_wavel/resolution/nPix_per_resol)

        # Call rotational and instrumental broadening mechanisms
        flux_rot = Apply_Rotation(wave_bin, flux_bin, vsini)
        flux_broad = Apply_Instrument(wave_bin, flux_rot, resolution)

        # Call rebinning for necessary parameters and get the total transmission
        mirror_reflectivity = Rebin_params(wave_color, mirror_ref_nobin, wave_bin, "quadratic")
        QE = Rebin_params(wave_color, QE_nobin, wave_bin, "quadratic")
        transmission = Get_transmission(wave_bin, mirror_reflectivity)

        # Blaze function of the spectrograph
        blaze_peak = Calc_Blaze_Peak(wave_bin, lambda_cen, lambda_min, lambda_max, order_numbers)

        # Call functions to extract observed simulated spectrum
        flux_eff = Calc_ConversionFlux(wave_bin, flux_broad, M1_area, QE, transmission, blaze_peak)
        flux, straylight = Calc_Flux_bin(wave_bin, flux_eff, ExpTime, resolution, transmission, QE)

        # Compute noise and SNR
        noise, snr = Calc_Noise(wave_bin, flux, straylight, ExpTime)

        SNR_ls.append(np.median(snr))

    plt.semilogy(mag_list, SNR_ls, label = "{} s".format(ExpTime))

# Stop the clock
end_time = datetime.now()

# Running time
print('Runtime:', end_time - start_time)

plt.axhline(200, ls = ':', color = 'k')
plt.xlim(mag_list[0], mag_list[-1])
plt.ylim(10, SNR_ls[0])
plt.xlabel('Magnitude')
plt.ylabel('SNR')
plt.legend(loc = 'best')
plt.grid(True, 'both')
plt.savefig(common_path+'/Figures/Teff{}g{}R{}D{}/SNR_MAG_T{}g{}D{}R{}vsini{}.pdf'.format(Teff, g, rad_sol, dist_pc, Teff, g, dist_pc, rad_sol, vsini), facecolor='w', transparent=False, format='pdf')
plt.show()