# Imports
import numpy as np

"""
Variables in code
"""

### Constants
cc = 299792458            # speed of light [m/s]
h = 6.626 * 1e-34         # Plank's constant [J*s]
sigma = 5.670 * 1e-8      # Stefan-Boltzmann constant [W/m²/K^4]
Rsun_pc = 2.25461 * 1e-8  # radius of the Sun [pc]

### Define the Teff, g, and v values you want to import (Tlusty BStarGrid2006)
Teff = 27000              # K (Beta Ceph)
g = 400                   # m/s² (Surface gravity)
v = 2                     # km/s (Microturbulent velocity)

### Parameters
rad_sol = 5.6             # R_sol (Beta Ceph)
dist_pc = 210             # pc
vsini = 100               # km/s (or 300)

### Stellar parameters
Vmag = 4                  # Magnitude
ExpTime = 900             # Exposure time (s)
SpecType = 'B'            # O, B, A, F, G, K, M

### V band zero point flux
flux_0_johnson_v = 3.735 * 1e-9 * 840 # erg/s/cm²

### Paths, filenames & directories
# Common path where the inputs and outputs are/will be stored
common_path = '/Users/philippen/Library/Mobile Documents/com~apple~CloudDocs/Desktop/PhD/LPV Simulator'
ster_path = '/STER/philippen/LPV_sim'

# Figures folder
fig_path = common_path + '/' + 'Figures/Teff{}g{}R{}D{}/'.format(Teff, g, rad_sol, dist_pc)

# Pulsation inputs
pulsation_path = common_path + '/' + 'Pulsations_famias'
# pulsation_path = ster_path + '/' + 'Pulsations_famias'

# Spectral orders .txt file
order_filename = 'blaze_orders.txt'

# Folder of the pulsation inputs
# pulsation_dir = 'beta_Cep_vsini100'
# pul_dirs = ['3months_15min_l0m0_11.53cd', '3months_15min_l1m0_11.53cd', '3months_15min_l1m1_11.53cd', '3months_15min_l8m0_11.53cd']
pul_dirs = ['3months_90min_l0m0_11.53cd', '3months_90min_l1m0_11.53cd', '3months_90min_l1m1_11.53cd', '3months_90min_l8m0_11.53cd']

# Folder of the static profile
null_profile_dir = 'beta_Cep_vsini100_static'

# Pulsation outputs (convolved spectra)
output_path = common_path + '/' + 'Pulsations_OUT/Teff{}g{}R{}D{}'.format(Teff, g, rad_sol, dist_pc)
# output_path = ster_path + '/' + 'Pulsations_OUT/Teff{}g{}R{}D{}'.format(Teff, g, rad_sol, dist_pc)

### CubeSPEC wavelength range (boundaries)
lambda_lower = 4200 # A
lambda_upper = 6200 # A
mean_wavel = 5200   # middle wavelength of CubeSPEC range [A]

### Spectral lines "data" lists of names, lower, upper and central wavelengths maybe?
# line_name = ['Hβ', 'SiIII', 'HeI']
# line_cen = [4861, 4552, 5876]
# line_inf = [4840, 4550, 5871]
# line_sup = [4882, 4555, 5881]
#
# line_name = ['HeI']
# line_cen = [5876]
# line_inf = [5871]
# line_sup = [5881]
#
# line_name = ['SiIII']
# line_cen = [4552]
# line_inf = [4550]
# line_sup = [4555]
#
line_name = ['SiIII', 'HeI']
line_cen = [4552, 5876]
line_inf = [4550, 5871]
line_sup = [4555, 5881]

noise_status = True

line_wvl = 5200 # For now the middle of the range for rebinning the pulsation profile for convolution matter (past)

### Spectrograph parameters
resolution = 55000              # R = lambda/delta_lambda of CubeSpec
grating_groove_density = 41.59  # lines/mm
grating_spacing = 24.0442414    # nm
blaze_angle = np.radians(76)    # radians
ratio_edge_centre_tm = 0.4

### CubeSPEC characteristics (wavelength independent)
# Intrinsic colors
UB = {'O': -1.20, 'B': -0.50, 'A': 0.10, 'F': 0.00, 'G': 0.20, 'K': 1.10, 'M': 1.30 }
UB = {'O': -0.80, 'B': -0.20, 'A': 0.10, 'F': 0.00, 'G': 0.10, 'K': 0.50, 'M': 0.70 }
BV = {'O': -0.30, 'B': -0.15, 'A': 0.20, 'F': 0.50, 'G': 0.70, 'K': 1.20, 'M': 1.70 }
VI = {'O': -0.15, 'B': -0.10, 'A': 0.15, 'F': 0.50, 'G': 0.76, 'K': 1.20, 'M': 3.00 }
RI = {'O': -0.10, 'B': -0.06, 'A': 0.08, 'F': 0.25, 'G': 0.36, 'K': 0.70, 'M': 1.80 }
ZI = {'O': -0.05, 'B': -0.03, 'A': 0.04, 'F': 0.10, 'G': 0.15, 'K': 0.40, 'M': 0.20 }

#Area = 19.2*8.3             # M1 area (cm²)
nMirrors = 9                 # Number of mirror reflections
disperser_eff = 0.67         # Efficiency of the main disperser (diffraction grating)
cross_disperser_eff = 0.90   # Efficiency of the cross disperser (1 = No cross dispersion)

M1_width = 7.5
M1_length = 18
Corner_area = 4 * 1.5
M1_area = M1_width * M1_length - Corner_area    # M1 area (cm²)

M2_width = 5.0
M2_length = 2.6
M2_area = M2_width * M2_length
obscuration = M2_area/M1_area               # Obscuration due to M2

FGS_reflectivity = 0.92                     # Flux loss for Fine Guide Sensor pick off or dichroic beam splitter
pointing_jitter_loss = 0.33                 #"""!!! # QUESTION !!!"""

read_noise = 2.0                            # Detector read noise per pixel (e-)
dark_flux = 1.3                             # Detector dark current per pixel (e-/s)
nPix_per_resol = 3.08                       # number of detector pixels per resol. elem. (slit is 20um wide and is sampled by 6.5um sized pixels: ~3.08)
nPixelsSpectral = 3                         # Spectral direction: 3 pixels/resolution element
nPixelsSpatial = 11                         # Spatial profile
nPixels = nPixelsSpectral * nPixelsSpatial  # Number of pixels (spectral x spatial direction)
gain_var = 0.0005                           # inter-pixel relative Gain Variability: gain is slightly different for each pixel

straylight_M1 = 170                         # photons/A/s (scattered from M1)
#straylight_M1 = 170 / nPixels              # photons/pixel/s (The straylight is "averaged" over the 33 pixels) %TODO: is that correct? NOPE
straylight_baffle = 1                       # photons/bin/s (reflected by baffle)
straylight_baffle = 1 / nPix_per_resol      # photons/pixel/s (3.08 pixels per resol. elem.)
straylight_internal = 0.03                  # fraction of local flux (Spectrograph straylight)

#### CubeSpec optical characteristics (wavelength dependent -> interpolation required)
# Colors wvl: ['U', 'B', 'V', 'R', 'I']
wave_color = [3800, 4400, 5500, 6400, 7900]  # Angstrom  #QUESTION: approximate, best use effective filter lam0?

# Mirror reflectivity in colors: ['U', 'B', 'V', 'R', 'I']
mirror_ref_nobin = [0.95, 0.97, 0.98, 0.98, 0.98]

# Quantum Efficiency in colors: ['U', 'B', 'V', 'R', 'I']
QE_nobin = [0.46, 0.72, 0.86, 0.84, 0.60]