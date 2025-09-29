import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import CubicSpline
import pandas as pd
import pickle
import glob
from scipy import signal
from barycorrpy import get_BC_vel , exposure_meter_BC_vel
from astropy.time import Time
from scipy.interpolate import interp1d
import plotly.express as px
from expecto import get_spectrum
from PyAstronomy import pyasl
import lmfit
import os
from lmfit import Model
import corner
import random
from datetime import datetime, timedelta
import parameters
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

base_path = '/Users/odin/Library/CloudStorage/OneDrive-UniversitaetBern/Code/RISTRETTO_Redux/5_MaddaelnaCodeBundle_250912/Josh/'

# Transmission function RISTRETTO, without Nico coupling maps 

T_0 = 96.6e-2 # Atmosphere transmission lost
T_1 = 61.2e-2 # Alluminium mirror coating
T_2 = 68.1e-2   # Front-end optical transmission
T_5 = 93.3e-2 # Raw Fiber Transmission
T_6 = 95.5e-2 # 3d printed Lens Transmission
T_7 = 43.9e-2 # Spectrograph efficiency

T_3 = 40e-2  # AO performances
T_4 = 40e-2 # Coronograph performances
T_34= 40e-2

partial_efficiency= T_0*T_1*T_2*T_5*T_6*T_7
print(f'Total transmission in optimal conditions: {T_0*T_1*T_2*T_5*T_6*T_7*100} %')

####################### IMPORT SYNTHETIC STELLAR SPECTRUM from BT-SETTL-CIFIST file (downsampled to R~140,000) ###########################
# flux in erg/s/cm2/A
    # Load your spectrum
wave_star_restFrame, spec_star_erg = np.loadtxt("/Users/odin/Library/CloudStorage/OneDrive-UniversitaetBern/Code/RISTRETTO_Redux/1_Get_BTSETTL_CIFIST_resample/PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt", unpack=True)
    # Create mask
mask = (wave_star_restFrame >= 6000) & (wave_star_restFrame <= 8600)
    # Apply mask to both arrays
wave_star_restFrame = wave_star_restFrame[mask]/10
spec_star_erg = spec_star_erg[mask]
    # convert flux to erg/s/cm2/cm
spec_starr = spec_star_erg * 1e8    

#### SAVE TO CSV
wavelength_str = [f"{w:.3f}" for w in wave_star_restFrame]
flux_str = [f"{f:.6e}" for f in spec_star_erg]  # scientific notation
data = np.column_stack((wavelength_str, flux_str))
np.savetxt("CSV_files/INPUT_btsettl_spectrum_PDS70_stars.csv", data, delimiter=",", header="wavelength_nm,flux", comments="", fmt="%s")
###

###
c                       = 2.99792458e8      # speed of light in m/s
h                       = 6.62607015e-34    # Planck's constant in J⋅Hz−1
pc                      = 3.08567758128e16  # 1 parsec, in m
sun_radius              = 695700000         # in m
stellar_radius          = 1.26*sun_radius   # in m 
stellar_distance        = 113               # in parsec
telescope_eff_surface   = 49.3              # Effective surface of the VLT in m^2, accounting for obstructions

# Transform the spectrum to photons/s/A seen on the VLT
spec_star             = spec_starr.copy()                                      # erg/s/cm2/cm
spec_star             = spec_star / 1e7                                        # J/s/cm2/cm
spec_star             = spec_star / (h * c / (wave_star_restFrame * 1e-9))     # photons/s/cm2/cm
spec_star             = spec_star * 1e4                                        # photons/s/m2/cm
spec_star             = spec_star * 4 * np.pi * stellar_radius**2              # photons/s/cm (total emitted flux, using PHOENIX radius, close to Proxima empirical luminosity of 5.93e30 erg/s)
spec_star             = spec_star / (4 * np.pi * (stellar_distance * pc)**2)   # photons/s/m2/cm (flux received at Earth)
spec_star             = spec_star * telescope_eff_surface                      # photons/s/cm  (flux entering the telescope)
spec_star             = spec_star*partial_efficiency/ 1e8                      # photons/s/angstrom
spec_star_ph          = spec_star.copy()

minim_wave=610
max_wavele=850

max_cut=150

sampling=0.35 #35      # km/s/pixel
logwave=np.arange(np.log10(minim_wave),np.log10(max_wavele),np.log10(1.+sampling/2.99792458e5))
wavelength_grid_angstroms=10**logwave
wave=wavelength_grid_angstroms[100:-max_cut]
spec_star_ph=(CubicSpline(wave_star_restFrame,spec_star_ph)(wavelength_grid_angstroms))

for fiber in [1]: #1, 2, 3, 4, 5, 6, 7, 8
    coupl_object=3.5e-4
    fiber_flux=CubicSpline(wavelength_grid_angstroms,spec_star_ph)(wavelength_grid_angstroms)*coupl_object
    f_fiber = f'fibertest_{fiber}.csv'
    csv_rows = ["{} , {}".format(wave[j],  fiber_flux[j]) for j in range(len(wave))]
    path_1= base_path + 'CSV_files/'
    with open( path_1+ f_fiber, 'w') as f:
        f.write("\n".join(csv_rows))
        
os.system(f'scp CSV_files/*.csv jb23l046@submit03.unibe.ch:/storage/homefs/jb23l046/Simu_run/data/in/')
