import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import textwrap
import smplotlib

# Ensure correct usage
if len(sys.argv) != 3:
    print("Usage: python add.py <file1> <file2>")
    sys.exit(1)

# Get filenames from command-line arguments
file1_name = sys.argv[1]
file2_name = sys.argv[2]

# Directory containing the input files
input_dir = "input_fromubelix/in"

# Build full file paths
file1_path = os.path.join(input_dir, file1_name + ".csv")
file2_path = os.path.join(input_dir, file2_name + ".csv")

# Load CSV files (assuming two columns: wavelength, flux)
wavelength1, flux1 = np.loadtxt(file1_path, delimiter=',', unpack=True)
wavelength2, flux2 = np.loadtxt(file2_path, delimiter=',', unpack=True)

# Check that wavelengths match (simple sanity check)
if not np.array_equal(wavelength1, wavelength2):
    print("Warning: Wavelength arrays do not match. Interpolation may be required.")
    sys.exit(1)

# Add the fluxes
flux_sum = flux1 + flux2

# Function to wrap long legend text
def wrap_label(label, width=40):
    return "\n".join(textwrap.wrap(label, width=width))

# Apply wrapping to legend labels
file1_label = wrap_label(file1_name)
file2_label = wrap_label(file2_name)
sum_label = wrap_label("Sum")

# Plot the three spectra
plt.figure(figsize=(10, 6))
plt.plot(wavelength1 * 10, flux1, label=file1_label, color='blue')
plt.plot(wavelength2 * 10, flux2, label=file2_label, color='green')
plt.plot(wavelength1 * 10, flux_sum, label=sum_label, color='red', linestyle='--')

plt.xlabel("Wavelength (Angstroms)")
plt.ylabel("Flux (photons/s)")
plt.title("INPUT: PDS70b 100mas separation off-axis external spaxel, incident at telescope")
plt.yscale('log')
plt.xlim(6550, 6575)
plt.grid(True)
plt.legend(loc='best', fontsize=9)  # keep it inside the plot
plt.tight_layout()
plt.show()
