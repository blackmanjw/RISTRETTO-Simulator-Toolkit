#!/usr/bin/env python3
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# --- Helper Functions ---

def gauss_convolve(spectrum, fwhm):
    sigma = fwhm / 2.3548
    offset = np.arange(-2*round(fwhm), 2*round(fwhm)+1)
    kernel = np.exp(-offset**2 / (2 * sigma**2))
    return np.convolve(spectrum, kernel / kernel.sum(), mode='same')

def create_spec_from_csv(df, exposure_time, wave, wave_resampled):
    resolution_rv = 0.86  # km/s
    flux = df['flux'].values * exposure_time * wave * 10 * resolution_rv / 2.99792458e5
    flux = gauss_convolve(flux, 300000 / 140000. / 0.35)
    return CubicSpline(wave, flux)(wave_resampled)

def add_noise(flux, RON2):
    noise = np.random.normal(0, np.sqrt(flux + RON2))
    return flux + noise

# --- Main script ---

if len(sys.argv) < 2:
    print("Usage: python script.py <filename.csv> [t_exp_seconds]")
    sys.exit(1)

filename = sys.argv[1]

# Optional t_exp argument
if len(sys.argv) >= 3:
    try:
        t_exp = float(sys.argv[2])
    except ValueError:
        print("t_exp must be a number (seconds).")
        sys.exit(1)
else:
    t_exp = 3600  # default 1 hour

# Automatically include the output folder
base_path = os.path.join(os.getcwd(), "output")
csv_file = os.path.join(base_path, filename)

if not os.path.exists(csv_file):
    print(f"File not found: {csv_file}")
    sys.exit(1)

# Load CSV with header
df_star = pd.read_csv(csv_file)

# Wavelength grids
minim_wave = 610.0
max_wavele = 850.0
sampling_raw = 0.35
sampling_resampled = 0.86
max_cut = 150

logwave = np.arange(np.log10(minim_wave), np.log10(max_wavele),
                    np.log10(1 + sampling_raw / 2.99792458e5))
wave = 10**logwave[100:-max_cut]

logwave2 = np.arange(np.log10(minim_wave), np.log10(max_wavele),
                     np.log10(1 + sampling_resampled / 2.99792458e5))
wave2 = 10**logwave2[100:-max_cut]

# Process spectrum
flux_star = create_spec_from_csv(df_star, t_exp, wave, wave2)

# Add noise
RON2 = 29
flux_star_noisy = add_noise(flux_star, RON2)

# Plot
plt.figure(figsize=(10,5))
plt.plot(wave2*10, flux_star_noisy, label=f'Noisy Spectrum (t_exp={t_exp}s)')
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (photons/s)')
plt.title(f"Noisy Spectrum from {filename}")
plt.legend()
plt.tight_layout()
plt.show()
