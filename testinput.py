import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from scipy.interpolate import CubicSpline

# Import Expecto for Phoenix spectrum
from expecto import get_spectrum  # adjust if your package name is different

# --------------------------
# Spectrum Loading Functions
# --------------------------
def load_file_spectrum(filename):
    """Load txt spectrum, ignoring header, and convert flux to cm units"""
    data = np.loadtxt(filename, comments='#')
    wavelength = data[:, 0]   # in Å
    flux = data[:, 1] * 1e8   # convert flux from erg/cm2/s/Å to erg/cm2/s/cm

    # Restrict to 6000-8500 Å
    mask = (wavelength >= 6000) & (wavelength <= 8500)
    wavelength = wavelength[mask] * 1e-1  # convert Å to nm
    flux = flux[mask]

    return wavelength, flux

def load_star_spectrum():
    """Generate Phoenix model spectrum in selected wavelength range"""
    spectrum_proxima = get_spectrum(T_eff=4100, log_g=5, cache=False, vacuum=True)
    
    val = np.where((spectrum_proxima.wavelength.value < 8600) &
                   (spectrum_proxima.wavelength.value > 6000))
    
    flux_star = spectrum_proxima.flux.value[val]
    wave_star = spectrum_proxima.wavelength.value[val] * 1e-1  # convert Å to nm
    return wave_star, flux_star

def plot_both_spectra(filename):
    """Plot observed and Phoenix spectra"""
    wave_file, flux_file = load_file_spectrum(filename)
    wave_star, flux_star = load_star_spectrum()

    plt.figure(figsize=(10,6))
    plt.plot(wave_file, flux_file, label=f'Observed: {filename}', color='blue', lw=1)
    plt.plot(wave_star, flux_star, label='Phoenix Model (4100K, logg=5)', color='red', lw=1)
    
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Flux (erg/cm²/s/cm)')
    plt.title('Observed vs Phoenix Star Spectrum')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# --------------------------
# Main Script
# --------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python spectrum_plot.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    #plot_both_spectra(filename)

    # Load observed spectrum
    wave_star_restFrame, spec_starr = load_file_spectrum(filename)
    
    # --------------------------
    # Transmission function RISTRETTO, without Nico coupling maps 
    # --------------------------
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


    # --------------------------
    # Physical constants
    # --------------------------
    c                       = 2.99792458e8      # m/s
    h                       = 6.62607015e-34    # J·s
    pc                      = 3.08567758128e16  # m
    sun_radius              = 695700000         # m
    stellar_radius          = 1.26*sun_radius   # m
    stellar_distance        = 113               # parsec
    telescope_eff_surface   = 49.3              # m^2

    # --------------------------
    # Transform spectrum to photons/s/Å
    # --------------------------
    spec_star             = spec_starr.copy()                                      # erg/s/cm2/cm
    spec_star             = spec_star / 1e7                                        # J/s/cm2/cm
    spec_star             = spec_star / (h * c / (wave_star_restFrame * 1e-9))     # photons/s/cm2/cm
    spec_star             = spec_star * 1e4                                        # photons/s/m2/cm
    spec_star             = spec_star * 4 * np.pi * stellar_radius**2              # photons/s/cm
    spec_star             = spec_star / (4 * np.pi * (stellar_distance * pc)**2)   # photons/s/m2/cm
    spec_star             = spec_star * telescope_eff_surface                      # photons/s/cm
    spec_star             = spec_star * partial_efficiency / 1e8                   # photons/s/Å
    spec_star_ph          = spec_star.copy()

    # --------------------------
    # Wavelength grid and resampling
    # --------------------------
    minim_wave = 610
    max_wavele = 850
    max_cut = 150
    sampling = 0.35  # km/s/pixel
    logwave = np.arange(np.log10(minim_wave), np.log10(max_wavele), np.log10(1. + sampling / 2.99792458e5))
    wavelength_grid_angstroms = 10**logwave
    wave = wavelength_grid_angstroms[100:-max_cut]
    spec_star_ph = CubicSpline(wave_star_restFrame, spec_star_ph)(wavelength_grid_angstroms)

    # Plot photon flux
    plt.figure(figsize=(10,5))
    plt.plot(wave_star_restFrame, spec_star_ph[:len(wave_star_restFrame)], color='red')
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Flux (photons/s/Å)")
    plt.title("Photon Flux after Telescope + Efficiency Corrections")
    plt.grid(True)
    plt.show()

    # --------------------------
    # Generate CSV files for fibers
    # --------------------------
    base_path = "./"  # adjust this to your desired folder
    for fiber in [1]:  # can extend to multiple fibers
        coupl_object = 3.5e-4
        fiber_flux = CubicSpline(wavelength_grid_angstroms, spec_star_ph)(wavelength_grid_angstroms) * coupl_object
        f_fiber = f'fiber_{fiber}.csv'
        csv_rows = ["{} , {}".format(wave[j], fiber_flux[j]) for j in range(len(wave))]
        path_1 = base_path + 'CSV_files/'
        import os
        os.makedirs(path_1, exist_ok=True)
        with open(path_1 + f_fiber, 'w') as f:
            f.write("\n".join(csv_rows))

    # --------------------------
    # Simulate noisy spectrum
    # --------------------------
    def gauss_convolve(spectrum, fwhm):
        sigma = fwhm / 2.3548
        offset = np.arange(-2*round(fwhm), 2*round(fwhm)+1)
        kernel = np.exp(-offset**2 / (2 * sigma**2))
        return np.convolve(spectrum, kernel / kernel.sum(), mode='same')

    def create_spec_from_csv(df, exposure_time, wave, wave_resampled):
        resolution_rv = 0.86  # km/s
        flux = df[1].values * exposure_time * wave * 10 * resolution_rv / 2.99792458e5
        flux = gauss_convolve(flux, 300000 / 140000. / 0.35)
        return CubicSpline(wave, flux)(wave_resampled)

    def add_noise(flux, RON2):
        noise = np.random.normal(0, np.sqrt(flux + RON2))
        return flux + noise

    # Wavelength grids
    sampling_raw = 0.35
    sampling_resampled = 0.86
    logwave2 = np.arange(np.log10(minim_wave), np.log10(max_wavele), np.log10(1 + sampling_resampled / 2.99792458e5))
    wave2 = 10**logwave2[100:-max_cut]

    # Load CSV and simulate
    t_exp = 3600  # 1 hour
    fiber_number = 1
    path_csv = base_path + 'CSV_files/'
    df_star = pd.read_csv(f'{path_csv}fiber_{fiber_number}.csv', header=None)
    flux_star = create_spec_from_csv(df_star, t_exp, wave, wave2)

    RON2 = 29
    flux_star_noisy = add_noise(flux_star, RON2)

    # Plot simulated noisy spectrum
    plt.figure(figsize=(10,5))
    plt.plot(wave2*10, flux_star_noisy, label='Noisy Spectrum')
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Flux (photons)')
    plt.title('Simulated Noisy Spectrum')
    plt.legend()
    plt.grid(True)
    plt.show()
