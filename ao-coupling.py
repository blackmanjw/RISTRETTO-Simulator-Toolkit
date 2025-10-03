#################
# AO Coupling Maps for RISTRETTO
# Features:
# - Average over phase-screen dimension
# - Smooth curves with moving average
# - Extra smooth curves using Gaussian smoothing

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

# --- Load FITS data ---
hdul = fits.open("ao-coupling/rho_offaxis_Imag_7.4_s_0.76_elev_90_wfsR_60_g_0.7_NmodeBR_0_f2_1.0_delay_1_Npix_120_misreg_0.00_0.25_-0.25_0.00_OG1.fits")
hdul.info()
data = hdul[0].data
distance = hdul[2].data.flatten()

# --- Helper functions ---
def avg_over_phase(x):
    """Average along phase-screen dimension if it exists."""
    return np.mean(x, axis=0) if x.ndim > 1 else x

def moving_average(y, window_size=7):
    """Smooth a 1D array with simple moving average."""
    return np.convolve(y, np.ones(window_size)/window_size, mode='same')

# --- Fibres ---
fibre1 = avg_over_phase(data[:, 2, 11, 0])
fibre2 = avg_over_phase(data[:, 2, 11, 1])
fibre3 = avg_over_phase(data[:, 2, 4, 2])
fibre4 = avg_over_phase(data[:, 2, 4, 3])
fibre5 = avg_over_phase(data[:, 2, 4, 4])
fibre6 = avg_over_phase(data[:, 2, 4, 5])
fibre7 = avg_over_phase(data[:, 2, 4, 6])

# Combine fibres 3–7 by averaging
fibre37 = (fibre3 + fibre4 + fibre5 + fibre6 + fibre7) / 5

# On-axis and off-axis fibres (for reference)
fibre_onaxis = avg_over_phase(data[:, 2, 2, 1])
fibre_offaxis = avg_over_phase(data[:, 2, 2, 0])

# Save fibre_offaxis vs distance to text
output = np.column_stack((distance, fibre_offaxis))
np.savetxt("ao-coupling/fibre_offaxis_vs_distance.txt", output,
           header="Distance  Fibre_OffAxis", fmt="%.8f")

# --- Initial smoothing with moving average ---
window_size = 7  # adjust for desired smoothness
fibre1_smooth = moving_average(fibre1, window_size)
fibre2_smooth = moving_average(fibre2, window_size)
fibre37_smooth = moving_average(fibre37, window_size)

# --- Save fibre2 smoothed vs distance to text ---
output_smoothed = np.column_stack((distance, fibre2_smooth))
np.savetxt("ao-coupling/fibre_offaxis_vs_distance.txt", output_smoothed,
           header="Distance  Fibre2_Smoothed", fmt="%.8f")
print("Saved ao-coupling/fibre_offaxis_vs_distance.txt")

# --- Extra smooth curves using Gaussian smoothing ---
sigma_extra = 5  # adjust for more/less extra smoothing
fibre1_extra = gaussian_filter1d(fibre1_smooth, sigma=sigma_extra)
fibre2_extra = gaussian_filter1d(fibre2_smooth, sigma=sigma_extra)
fibre37_extra = gaussian_filter1d(fibre37_smooth, sigma=sigma_extra)

# --- Plotting ---
plt.figure(figsize=(8,5))

# Smoothed curves
plt.plot(distance, fibre1_smooth, label="Fibre 1", color='#0072bd')
plt.plot(distance, fibre2_smooth, label="Fibre 2", color='#d95319')
plt.plot(distance, fibre37_smooth, label="Fibre 3-7 (avg)", color='#edb120')

# Extra smooth curves (dashed)
#plt.plot(distance, fibre1_extra, color='#0072bd', linestyle='--', linewidth=2)
#plt.plot(distance, fibre2_extra, color='#d95319', linestyle='--', linewidth=2)
#plt.plot(distance, fibre37_extra, color='#edb120', linestyle='--', linewidth=2)

plt.xlabel("Off-axis distance (arcseconds)")
plt.ylabel("Stellar Coupling")
plt.yscale("log")
plt.xlim(-0.02, 0.7)
plt.legend()
plt.tight_layout()
plt.savefig('ao-coupling/ao-coupling.png', dpi=300)
plt.show()
