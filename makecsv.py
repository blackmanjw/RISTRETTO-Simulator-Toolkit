import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, interp1d
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

############################
# Handle arguments
############################
if len(sys.argv) < 2 or len(sys.argv) > 3:
    print("Usage: python makecsv.py filename.txt [separation_mas]")
    sys.exit(1)

filename = sys.argv[1]
separation_mas = None
if len(sys.argv) == 3:
    separation_mas = float(sys.argv[2])
    separation_arcsec = separation_mas / 1000.0  # convert mas to arcsec

# Input directory
input_dir = os.path.join(os.getcwd(), "input")
input_file = os.path.join(input_dir, filename)

if not os.path.exists(input_file):
    print(f"Error: {input_file} not found.")
    sys.exit(1)

# Strip extension for naming
filename_base = os.path.splitext(os.path.basename(filename))[0]

# Make sure output dir exists
output_dir = os.path.join(os.getcwd(), "output")
os.makedirs(output_dir, exist_ok=True)

############################
# Constants
############################
T_0 = 96.6e-2
T_1 = 61.2e-2
T_2 = 68.1e-2
T_5 = 93.3e-2
T_6 = 95.5e-2
T_7 = 43.9e-2

partial_efficiency = T_0 * T_1 * T_2 * T_5 * T_6 * T_7
print(f"Total transmission in optimal conditions: {partial_efficiency*100:.2f} %")

############################
# Load spectrum from file
############################
wave_star_restFrame, spec_star_erg = np.loadtxt(input_file, unpack=True)

mask = (wave_star_restFrame >= 6000) & (wave_star_restFrame <= 8600)
wave_star_restFrame = wave_star_restFrame[mask] / 10
spec_star_erg = spec_star_erg[mask]

spec_starr = spec_star_erg * 1e8

############################
# Save to CSV (raw spectrum)
############################
wavelength_str = [f"{w:.3f}" for w in wave_star_restFrame]
flux_str = [f"{f:.6e}" for f in spec_star_erg]
data = np.column_stack((wavelength_str, flux_str))

csv_outfile = os.path.join(output_dir, f"{filename_base}.csv")
np.savetxt(csv_outfile, data, delimiter=",", header="wavelength_nm,flux", comments="", fmt="%s")

############################
# Photon flux calculation
############################
c = 2.99792458e8
h = 6.62607015e-34
pc = 3.08567758128e16
sun_radius = 695700000
stellar_radius = 1.26 * sun_radius
stellar_distance = 113
telescope_eff_surface = 49.3

spec_star = spec_starr.copy()
spec_star = spec_star / 1e7
spec_star = spec_star / (h * c / (wave_star_restFrame * 1e-9))
spec_star = spec_star * 1e4
spec_star = spec_star * 4 * np.pi * stellar_radius**2
spec_star = spec_star / (4 * np.pi * (stellar_distance * pc)**2)
spec_star = spec_star * telescope_eff_surface
spec_star = spec_star * partial_efficiency / 1e8
spec_star_ph = spec_star.copy()

############################
# Resampling
############################
minim_wave = 610
max_wavele = 850
max_cut = 150
sampling = 0.35

logwave = np.arange(
    np.log10(minim_wave),
    np.log10(max_wavele),
    np.log10(1.0 + sampling / 2.99792458e5),
)
wavelength_grid_angstroms = 10**logwave
wave = wavelength_grid_angstroms[100:-max_cut]

spec_star_ph = PchipInterpolator(wave_star_restFrame, spec_star_ph)(wavelength_grid_angstroms)

############################
# Determine coupling object
############################
# Default values
coup_object = 0.45  # default
append_str = ""

if "star" in filename.lower():
    coup_object = 3.5e-4  # default for star
    print("This is a star.")

    if separation_mas is not None:
        separation_arcsec = separation_mas / 1000.0
        offaxis_file = "ao-coupling/fibre_offaxis_vs_distance_pchip.txt"
        print(f"This is an off-axis observation with separation: {separation_mas:.0f} mas")
        if os.path.exists(offaxis_file):
            offaxis_data = np.loadtxt(offaxis_file)
            distances = offaxis_data[:, 0]
            couplings = offaxis_data[:, 1]
            interp_func = interp1d(distances, couplings, bounds_error=False, fill_value="extrapolate")
            coup_object = float(interp_func(separation_arcsec))
        separation_str = f"{int(separation_mas)}mas"
    else:
        print("This is an on-axis observation.")
        separation_str = "38mas"

else:
    print("This is a planet.")
    separation_str = ""

# Construct append_str: coup_object first, then separation
if 0.4 <= coup_object <= 0.6:
    print(f"Adaptive Optics Coupling: {coup_object:.2f} %")
    coup_str = f"{coup_object:.2f}"
else:
    print(f"Adaptive Optics Coupling: {coup_object:.3e} %")
    coup_str = f"{coup_object:.2e}"

# Combine coup_object and separation
append_str = f"{coup_str}_{separation_str}" if separation_str else coup_str



print(f"Saved raw spectrum CSV -> {csv_outfile}")

############################
# Per-fiber output + plot
############################
for fiber in [1]:
    fiber_flux = (
        CubicSpline(wavelength_grid_angstroms, spec_star_ph)(wavelength_grid_angstroms)
        * coup_object
    )

    fiber_outfile = os.path.join(output_dir, f"{filename_base}_fiber{fiber}_{append_str}.csv")
    with open(fiber_outfile, "w") as f:
        f.write("wavelength_nm,flux\n")
        for j in range(len(wave)):
            f.write(f"{wave[j]:.3f},{fiber_flux[j]:.6e}\n")

    print(f"Saved fiber output -> {fiber_outfile}")

    # Save plot as PNG
    plt.figure(figsize=(8,5))
    plt.plot(wave, fiber_flux[:len(wave)], label=f"Fiber {fiber}")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Flux (photons/s)")
    plt.yscale("log")
    plt.title(f"Spectrum for {filename_base} - Fiber {fiber}")
    plt.legend()
    plt.tight_layout()

    png_outfile = os.path.join(output_dir, f"{filename_base}_fiber{fiber}_{append_str}.png")
    plt.savefig(png_outfile, dpi=150)
    plt.close()
    print(f"Saved plot -> {png_outfile}")
