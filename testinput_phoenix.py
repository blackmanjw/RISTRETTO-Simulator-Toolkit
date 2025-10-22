# %%
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
from lmfit import Model
import corner
import random
from datetime import datetime, timedelta
import parameters
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

# %%
base_path = '/Users/odin/Library/CloudStorage/OneDrive-UniversitaetBern/Code/RISTRETTO_Redux/5_MaddaelnaCodeBundle_250912/Josh/'

# %%
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

# %%
# Expecto package to create your star spectrum (Phoenix model)
spectrum_proxima = get_spectrum(T_eff=4100, log_g=5, cache=False, vacuum=True)
val = np.where((spectrum_proxima.wavelength.value < 8600) & (spectrum_proxima.wavelength.value > 6000))
spec_starr = spectrum_proxima.flux.value[val]
wave_star_restFrame = spectrum_proxima.wavelength.value[val] * 1e-1  # in nm

# --- Plot 1: Original PHOENIX spectrum ---
plt.figure(figsize=(10,5))
plt.plot(wave_star_restFrame, spec_starr, color='blue')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Flux (erg/s/cm²/cm)")
plt.title("Original PHOENIX Spectrum")
plt.grid(True)
plt.show()

# %%
# Physical constants
c                       = 2.99792458e8      # speed of light in m/s
h                       = 6.62607015e-34    # Planck's constant in J⋅Hz−1
pc                      = 3.08567758128e16  # 1 parsec, in m
sun_radius              = 695700000         # in m
stellar_radius          = 1.26*sun_radius   # in m 
stellar_distance        = 113               # in parsec
telescope_eff_surface   = 49.3              # Effective surface of the VLT in m^2, accounting for obstructions

# %%
# Transform the spectrum to photons/s/A seen on the VLT
spec_star             = spec_starr.copy()                                      # erg/s/cm2/cm
spec_star             = spec_star / 1e7                                        # J/s/cm2/cm
spec_star             = spec_star / (h * c / (wave_star_restFrame * 1e-9))     # photons/s/cm2/cm
spec_star             = spec_star * 1e4                                        # photons/s/m2/cm
spec_star             = spec_star * 4 * np.pi * stellar_radius**2              # photons/s/cm (total emitted flux)
spec_star             = spec_star / (4 * np.pi * (stellar_distance * pc)**2)   # photons/s/m2/cm (flux received at Earth)
spec_star             = spec_star * telescope_eff_surface                      # photons/s/cm  (flux entering the telescope)
spec_star             = spec_star * partial_efficiency / 1e8                   # photons/s/angstrom
spec_star_ph          = spec_star.copy()

minim_wave = 610
max_wavele = 850
max_cut = 150
sampling = 0.35  # km/s/pixel
logwave = np.arange(np.log10(minim_wave), np.log10(max_wavele), np.log10(1. + sampling / 2.99792458e5))
wavelength_grid_angstroms = 10**logwave
wave = wavelength_grid_angstroms[100:-max_cut]
spec_star_ph = CubicSpline(wave_star_restFrame, spec_star_ph)(wavelength_grid_angstroms)

# --- Plot 2: Photon flux after telescope + efficiencies ---
plt.figure(figsize=(10,5))
plt.plot(wave_star_restFrame, spec_star_ph[:len(wave_star_restFrame)], color='red')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Flux (photons/s/Å)")
plt.title("Photon Flux after Telescope + Efficiency Corrections")
plt.grid(True)
plt.show()

# %%
for fiber in [1]:  # can extend to 2..8
    coupl_object = 3.5e-4
    fiber_flux = CubicSpline(wavelength_grid_angstroms, spec_star_ph)(wavelength_grid_angstroms) * coupl_object
    f_fiber = f'fiber_{fiber}.csv'
    csv_rows = ["{} , {}".format(wave[j],  fiber_flux[j]) for j in range(len(wave))]
    path_1 = base_path + 'CSV_files/'
    with open(path_1 + f_fiber, 'w') as f:
        f.write("\n".join(csv_rows))

# %% [markdown]
# # Test detection flux ratio

# %%
waveSol = pickle.load(open(base_path + 'pickleFilesMaddalena/waveSol_from_HDF.pickle', "rb"))

def get_wave(waveSol, the_fiber, the_order):
    wavelength = waveSol[(waveSol['fiber_ID'] == the_fiber) & (waveSol['order_ID'] == the_order)]['poly'].values[0](np.arange(4096))
    return wavelength * 1e4

flat_pickle = pickle.load(open('pickleFilesMaddalena/master5.pickle', "rb"))

# %%
# --- Helper Functions ---

def gauss_convolve(spectrum, fwhm):
    """Convolve spectrum with a Gaussian kernel of given FWHM (in pixels)."""
    sigma = fwhm / 2.3548
    offset = np.arange(-2*round(fwhm), 2*round(fwhm)+1)
    kernel = np.exp(-offset**2 / (2 * sigma**2))
    return np.convolve(spectrum, kernel / kernel.sum(), mode='same')

def create_spec_from_csv(df, exposure_time, wave, wave_resampled):
    """Create simulated spectrum from CSV flux, convolved and interpolated."""
    resolution_rv = 0.86  # km/s
    flux = df[1].values * exposure_time * wave * 10 * resolution_rv / 2.99792458e5
    flux = gauss_convolve(flux, 300000 / 140000. / 0.35)
    return CubicSpline(wave, flux)(wave_resampled)

def add_noise(flux, RON2):
    """Add Poisson + readout noise."""
    noise = np.random.normal(0, np.sqrt(flux + RON2))
    return flux + noise

# --- Setup Wavelength Grids ---

sampling_raw = 0.35
sampling_resampled = 0.86
max_cut = 150
logwave = np.arange(np.log10(minim_wave), np.log10(max_wavele), np.log10(1 + sampling_raw / 2.99792458e5))
wave = 10**logwave[100:-max_cut]

logwave2 = np.arange(np.log10(minim_wave), np.log10(max_wavele), np.log10(1 + sampling_resampled / 2.99792458e5))
wave2 = 10**logwave2[100:-max_cut]

# --- Load and Process CSV Flux ---

t_exp = 1 * 3600        # Exposure time in seconds
fiber_number = 1

path_csv = base_path + 'CSV_files/'
df_star = pd.read_csv(f'{path_csv}fiber_{fiber_number}.csv', header=None)

flux_star = create_spec_from_csv(df_star, t_exp, wave, wave2)

# --- Add noise and plot ---

RON2 = 29  # Readout noise^2
flux_star_noisy = add_noise(flux_star, RON2)

plt.figure(figsize=(10,5))
plt.plot(wave2*10, flux_star_noisy, label='Noisy Spectrum')
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (photons/s)')
plt.title('Simulated Noisy Spectrum')
plt.legend()
plt.grid(True)
plt.show()

