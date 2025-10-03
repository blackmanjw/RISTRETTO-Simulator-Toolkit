## THIS uses the digitized AO coupling map from N. Billot and smoothes it for input for readcsv.txt The output is fibre_offaxis_vs_distance_lowess.txt.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import PchipInterpolator

# --------------------
# Read CSV
# --------------------
file_path = "ao-coupling/fibre2_digitized.csv"
data = pd.read_csv(file_path, header=None, names=['wavelength', 'star_coupling'])
data = data.sort_values(by='wavelength')

# Avoid zero for log-scale
epsilon = 1e-6
x = data['wavelength'].values
y = data['star_coupling'].values + epsilon

# --------------------
# LOWESS smoothing
# --------------------
frac = 0.1  # Adjust for smoothness
lowess_smoothed = lowess(y, x, frac=frac)
x_lowess, y_lowess = lowess_smoothed[:, 0], lowess_smoothed[:, 1]

# --------------------
# Save LOWESS fit
# --------------------
output_file = "ao-coupling/fibre_offaxis_vs_distance_lowess.txt"
np.savetxt(output_file, np.column_stack([x_lowess, y_lowess]),
           header="# Distance  Fibre_OffAxis", comments='', fmt='%.8e')
print(f"LOWESS fit saved to {output_file}")

# --------------------
# PCHIP interpolation
# --------------------
pchip = PchipInterpolator(x, y)
x_pchip = np.linspace(x.min(), x.max(), 1000)
y_pchip = pchip(x_pchip)

# --------------------
# Plot
# --------------------
plt.figure(figsize=(9,5))
plt.plot(x, y, 'o', label='Data', markersize=4)
plt.plot(x_lowess, y_lowess, '-', label=f'LOWESS Trend (frac={frac})', linewidth=2)
plt.plot(x_pchip, y_pchip, '--', label='PCHIP Interpolator', linewidth=2)
plt.yscale('log')
plt.xlabel("Wavelength")
plt.ylabel("Star Coupling")
plt.title("Star Coupling vs Wavelength: LOWESS vs PCHIP")
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()
