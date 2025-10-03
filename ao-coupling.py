#################
# THIS CODE IS TO PLOT THE AO COUPLING MAPS FOR RISTRETTO 
# using a .fits file supplied by Nico Blind.

from astropy.io import fits
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d 
import smplotlib

hdul = fits.open("ao-coupling/rho_offaxis_Imag_7.4_s_0.76_elev_90_wfsR_60_g_0.7_NmodeBR_0_f2_1.0_delay_1_Npix_120_misreg_0.00_0.25_-0.25_0.00_OG1.fits")
hdul.info()
data = hdul[0].data
distance = hdul[2].data
distance = distance.flatten()
fibre1 = data[:, 2, 11, 0]
fibre2 = data[:, 2, 11, 1]
fibre3 = data[:, 2, 4, 2]
fibre4 = data[:, 2, 4, 3]
fibre5 = data[:, 2, 4, 4]
fibre6 = data[:, 2, 4, 5]
fibre7 = data[:, 2, 4, 6]
fibre37=(fibre3+fibre4+fibre5+fibre6+fibre7)/5

# Apply Gaussian smoothing (adjust sigma for smoothness)
fibre1_smooth = gaussian_filter1d(fibre1, sigma=1.6)
fibre2_smooth = gaussian_filter1d(fibre2, sigma=1.6)
fibre37_smooth = gaussian_filter1d(fibre37, sigma=2)

fibre_onaxis = data[:, 2, 2, 1]
fibre_offaxis = data[:, 2, 2, 0]
distance[0:0]
#print(fibre_onaxis, distance)

# Stack distance and fibre_onaxis into 2 columns
output = np.column_stack((distance, fibre_offaxis))

# Save to text file
np.savetxt("fibre_offaxis_vs_distance.txt", output, header="Distance  Fibre_OffAxis", fmt="%.8f" )

# --- Interpolation function ---
def get_fibre_values_at(d):
    """Return interpolated fibre_onaxis and fibre_offaxis at distance d."""
    f_on = np.interp(d, distance, fibre_onaxis)
    f_off = np.interp(d, distance, fibre_offaxis)
    return f_on, f_off

# Example: query at distance = 0.5
query_distance = 0.1768
on_val, off_val = get_fibre_values_at(query_distance)
print(f"At distance {query_distance:.4f}: On-axis={on_val:.8e}, Off-axis={off_val:.8e}")

cmap = plt.get_cmap("viridis") 

plt.plot(distance, fibre1, label="Fibre 1", color='#440154FF')
plt.plot(distance, fibre2, label="Fibre 2", color='#2A788EFF')
#plt.plot(distance, fibre37, label="Fibre 3-7 (avg)", color='y')
plt.plot(distance, fibre37, label="Fibre 3-7 (avg)", color='#7AD151FF')
plt.xlabel("Off-axis distance (arcseconds)")
plt.ylabel("Stellar Coupling")
plt.yscale("log")
plt.xlim(-0.02,0.7)
plt.legend()
plt.savefig('coupling.png', bbox_inches='tight', dpi=300)
plt.show()
