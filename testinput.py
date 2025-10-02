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

t_exp = 1*3600        # Exposure time in seconds
fiber_number = 1

path_csv = base_path + 'output/'
df_star = pd.read_csv(f'{path_csv}fiber_{fiber_number}.csv', header=None)

flux_star = create_spec_from_csv(df_star, t_exp, wave, wave2)

# --- Add noise and loop over orders ---

RON2 = 29  # Readout noise^2
flux_star_noisy = add_noise(flux_star, RON2)

plt.plot(wave2*10, flux_star_noisy, label='Noisy Spectrum')
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (photons/s)')
plt.show()
