import numpy as np
import matplotlib.pyplot as plt
import os

# File path
input_file = "ao-coupling/fibre_offaxis_vs_distance.txt"

# Load data, skip header starting with #
data = np.loadtxt(input_file, comments="#")
distance = data[:, 0]  # in arcseconds
coupling = data[:, 1]  # fibre off-axis coupling

# Plot
plt.figure(figsize=(8,5))
plt.plot(distance, coupling, marker='o', linestyle='-')
plt.xlabel("Distance (arcsec)")
plt.ylabel("Fibre Off-Axis Coupling")
plt.title("Fibre Coupling vs Off-Axis Distance")
plt.grid(True)
plt.tight_layout()
plt.yscale("log")

# Output file in the same folder, same base name as PNG
output_file = os.path.splitext(input_file)[0] + ".png"
plt.savefig(output_file, dpi=150)
plt.close()

print(f"Saved plot -> {output_file}")
