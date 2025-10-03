#!/usr/bin/env python3
import subprocess
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import re
import smplotlib

# --- Handle command-line arguments ---
if len(sys.argv) != 2:
    print("Usage: python run_makecsv_and_plot.py <input_file.txt>")
    sys.exit(1)

input_arg = sys.argv[1]

# If input_arg has no directory, prepend input/
if os.path.dirname(input_arg) == "":
    input_file = os.path.join("input", input_arg)
else:
    input_file = input_arg

if not os.path.exists(input_file):
    print(f"Error: {input_file} not found.")
    sys.exit(1)

output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Extract base name to filter CSV files later
base_name = os.path.splitext(os.path.basename(input_file))[0]

# --- Step 1: Run makecsv.py for separations 100-700 mas in 50 mas steps ---
separations = range(100, 650, 50)
for sep in separations:
    print(f"Running makecsv.py with separation {sep} mas...")
    
    # Use the relative path for makecsv.py: just the filename in input/
    subprocess.run(["python", "makecsv.py", os.path.basename(input_file), str(sep)], check=True)


# --- Step 2: Collect generated CSV files for this base_name ---
csv_files = []
for f in os.listdir(output_dir):
    if f.startswith(base_name) and f.endswith(".csv") and re.search(r'(\d+)mas', f):
        csv_files.append(os.path.join(output_dir, f))

# Sort files by separation
csv_files.sort(key=lambda x: int(re.search(r'(\d+)mas', x).group(1)))

# --- Step 3: Plot all CSVs on the same axis with inline labels ---
plt.figure(figsize=(10,6))

# Keep track of previous label positions to reduce overlap
y_positions = []

for csv_file in csv_files:
    data = np.loadtxt(csv_file, delimiter=",", skiprows=1)
    wavelength = data[:, 0]
    flux = data[:, 1]

    # Extract separation label
    sep_label = re.search(r'(\d+)mas', csv_file).group(1) + " mas"

    # Plot line
    plt.plot(wavelength, flux, linewidth=1.5)

    # Determine inline label position (avoid overlapping previous labels)
    x_text = wavelength[-50]  # slightly before the end
    y_text = flux[-50]
    offset = 0.5
    while any(abs(y_text - yp) / yp < 0.05 for yp in y_positions):
        y_text *= (1 + offset)
        offset *= 1.1
    y_positions.append(y_text)

    plt.text(x_text, y_text, sep_label, fontsize=9, verticalalignment='bottom')

plt.xlabel("Wavelength (nm)")
plt.ylabel("Flux (photons/s)")
plt.title(f"{base_name} Spectra at Different Off-Axis Separations")
plt.yscale("log")
plt.tight_layout()

# Save figure
plot_file = os.path.join(output_dir, f"{base_name}_all_separations.png")
plt.savefig(plot_file, dpi=150)
plt.show()
print(f"Saved combined plot -> {plot_file}")
