import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import smplotlib
import os

path = "/Users/odin/Library/CloudStorage/OneDrive-UniversitaetBern/Code/RISTRETTO_Redux/1_Get_BTSETTL_CIFIST_resample/"   # <-- change to your desired path
filename = "PDS70b_BT-Settl-CIFIST-1400K-4logg.txt"

# Load your spectrum
wavelength, flux = np.loadtxt(os.path.join(path, filename), unpack=True)

# Desired resolution
R_target = 140000  # change as needed

# Compute delta_lambda at each point (approximate)
delta_lambda = wavelength / R_target
sigma_pixels = delta_lambda / (wavelength[1] - wavelength[0]) / 2.355

# If the wavelength grid is uniform, you can take the median value
sigma_median = np.median(sigma_pixels)

# Convolve flux with Gaussian kernel
flux_smoothed = gaussian_filter1d(flux, sigma=sigma_median)

# Plot result
plt.plot(wavelength, flux, label="Original")
plt.plot(wavelength, flux_smoothed, label=f"Downsampled (R={R_target})")
plt.legend()
plt.xlim(6200,8400)
#plt.ylim(0.00001,1000)
plt.ylim(0.00001,3e3)
plt.yscale('log')
plt.xlabel("Wavelength")
plt.ylabel("Flux")
# Save to file
name, ext = os.path.splitext(filename)
filename_save = f"{name}_140000{ext}"
plt.savefig(os.path.join("output", f"{type}{output_name}.png"), bbox_inches='tight', dpi=300)
plt.show()

np.savetxt(os.path.join(path, filename_save), np.column_stack([wavelength, flux_smoothed]))
