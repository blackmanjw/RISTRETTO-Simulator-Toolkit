import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator, interp1d
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import smplotlib

############################
# Handle arguments
############################
if len(sys.argv) < 2 or len(sys.argv) > 3:
    print("Usage: python makecsv.py filename.txt [separation_mas]")
    sys.exit(1)

filename = sys.argv[1]
separation_mas = None

# Parse optional argument (separation)
if len(sys.argv) == 3:
    try:
        separation_mas = float(sys.argv[2])
        separation_arcsec = separation_mas / 1000.0
    except ValueError:
        print(f"Warning: Unrecognized argument '{sys.argv[2]}' (ignored).")

# Input directory
input_dir = os.path.join(os.getcwd(), "input")
input_file = os.path.join(input_dir, filename)
if not os.path.exists(input_file):
    print(f"Error: {input_file} not found.")
    sys.exit(1)

filename_base = os.path.splitext(os.path.basename(filename))[0]

# Output directories
output_dir = os.path.join(os.getcwd(), "output")
os.makedirs(output_dir, exist_ok=True)
raw_dir = os.path.join(output_dir, "raw")
os.makedirs(raw_dir, exist_ok=True)

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
print(f"Total transmission in optimal conditions: {partial_efficiency*100:.2f} %")

############################
# Load spectrum
############################
wave_star_restFrame, spec_star_erg = np.loadtxt(input_file, unpack=True)
mask = (wave_star_restFrame >= 6000) & (wave_star_restFrame <= 8600)
wave_star_restFrame = wave_star_restFrame[mask] / 10  ## convert from A to nm
spec_star_erg = spec_star_erg[mask] # erg/s/cm2/A
spec_star_erg = spec_star_erg * 1e8 # erg/s/cm2/cm

############################
# Save raw spectrum
############################
wavelength_str = [f"{w:.3f}" for w in wave_star_restFrame]
flux_str = [f"{f:.6e}" for f in spec_star_erg]
data = np.column_stack((wavelength_str, flux_str))

csv_outfile = os.path.join(raw_dir, f"{filename_base}.csv")
np.savetxt(csv_outfile, data, delimiter=",", header="wavelength_nm,flux", comments="", fmt="%s")

############################
# Photon flux calc
############################
c = 2.99792458e8
h = 6.62607015e-34
pc = 3.08567758128e16
sun_radius = 695700000
stellar_distance = 113
telescope_eff_surface = 49.3
sun_radius = 695700000             # meters
jupiter_radius = 69911e3           # meters

########## SET RADIUS ACCORDING TO INPUT FILE

# Default radius (PDS70 star)
stellar_radius = 1.26 * sun_radius
radius_label = "1.26 R☉ (default)"

# Determine object type from filename
fname_lower = filename.lower()

if "pds70_star" in fname_lower:
    stellar_radius = 1.26 * sun_radius
    radius_label = "1.26 R☉"
elif "pds70b" in fname_lower:
    stellar_radius = 2.0 * jupiter_radius
    radius_label = "2.0 R_Jup"
elif "pds70c" in fname_lower:
    stellar_radius = 1.6 * jupiter_radius
    radius_label = "1.6 R_Jup"
elif "wispit2b" in fname_lower:
    stellar_radius = 1.6 * jupiter_radius
    radius_label = "1.6 R_Jup"
elif "2mj1612b" in fname_lower:
    stellar_radius = 1.5 * jupiter_radius
    radius_label = "1.5 R_Jup"
elif "2mj1612_star" in fname_lower:
    stellar_radius = 1.2 * sun_radius
    radius_label = "1.2 R☉"
elif "wispit2_star" in fname_lower:
    stellar_radius = 1.418 * sun_radius
    radius_label = "1.418 R☉"

print(f"→ Stellar radius set to {stellar_radius:.3e} m ({radius_label}) for '{filename_base}'")

###############################

spec_star             = spec_star_erg / 1e7                                    # J/s/cm2/cm
spec_star             = spec_star / (h * c / (wave_star_restFrame * 1e-9))     # photons/s/cm2/cm
spec_star             = spec_star * 1e4                                        # photons/s/m2/cm
spec_star             = spec_star * 4 * np.pi * stellar_radius**2              # photons/s/cm (total emitted flux, using PHOENIX radius, close to Proxima empirical luminosity of 5.93e30 erg/s)
spec_star             = spec_star / (4 * np.pi * (stellar_distance * pc)**2)   # photons/s/m2/cm (flux received at Earth)
spec_star             = spec_star * telescope_eff_surface                      # photons/s/cm  (flux entering the telescope)
spec_star             = spec_star*partial_efficiency/ 1e8                      # photons/s/angstrom
spec_star_ph          = spec_star.copy()

# Plot photon flux
plt.figure(figsize=(10,5))
plt.plot(wave_star_restFrame, spec_star_ph[:len(wave_star_restFrame)], color='red')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Flux (photons/s/Å)")
plt.title("Photon Flux after Telescope + Efficiency Corrections")
plt.grid(True)
plt.ticklabel_format(style='plain')
plt.savefig("test.png", dpi=150)

############################
# Plot flux at telescope (erg/s/cm²/Å)
############################
# Convert from stellar surface flux (spec_star_erg) to flux at telescope
# Apply scaling by (R*/d)²
flux_at_telescope = spec_star_erg * (stellar_radius / (stellar_distance * pc))**2

# Plot and save
plt.figure(figsize=(8, 5))
plt.plot(wave_star_restFrame, flux_at_telescope, color='tab:blue', lw=1)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Flux (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)")
plt.title(f"Flux at Telescope for {filename_base}")
plt.ticklabel_format(style='plain')
plt.yscale("log")
plt.tight_layout()

flux_plot_outfile = os.path.join(input_dir, f"{filename_base}_fluxattelescope.png")
plt.savefig(flux_plot_outfile, dpi=150)
plt.close()

print(f"Saved flux-at-telescope plot -> {flux_plot_outfile}")

############################
# Resampling
############################
minim_wave, max_wavele, max_cut, sampling = 610, 850, 150, 0.35
logwave = np.arange(np.log10(minim_wave), np.log10(max_wavele),
                    np.log10(1.0 + sampling / 2.99792458e5))
wavelength_grid_angstroms = 10**logwave
wave = wavelength_grid_angstroms[100:-max_cut]
spec_star_ph = PchipInterpolator(wave_star_restFrame, spec_star_ph)(wavelength_grid_angstroms)

############################
# Coupling logic
############################
# Default coupling list: always include 0.45
coupling_list = [0.45]

# If it's a star and no separation, also run 3.5e-4
if "star" in filename.lower() and separation_mas is None:
    coupling_list.append(3.5e-4)

# If separation specified, compute coupling from file (overrides others)
if separation_mas is not None:
    offaxis_file = "ao-coupling/fibre_offaxis_vs_distance_pchip.txt"
    if os.path.exists(offaxis_file):
        offaxis_data = np.loadtxt(offaxis_file)
        distances = offaxis_data[:, 0]
        couplings = offaxis_data[:, 1]
        interp_func = interp1d(distances, couplings, bounds_error=False, fill_value="extrapolate")
        coup_value = float(interp_func(separation_mas / 1000.0))
        coupling_list = [coup_value]
        print(f"Off-axis observation at {separation_mas:.0f} mas -> coupling {coup_value:.3e}")
    else:
        print("Warning: AO coupling file not found. Using default 3.5e-4")
        coupling_list = [3.5e-4]

############################
# Run for each coupling
############################
for coup_object in coupling_list:
    if 0.4 <= coup_object <= 0.6:
        coup_str = f"{coup_object:.2f}"
    else:
        coup_str = f"{coup_object:.2e}"

    append_str = f"{coup_str}"
    if separation_mas is not None:
        append_str = f"{coup_str}_{int(separation_mas)}mas"

    fiber_outfile = os.path.join(output_dir, f"{filename_base}_{append_str}.csv")
    png_outfile = os.path.join(output_dir, f"{filename_base}_{append_str}.png")

    print("\n============================")
    print(f"Coupling: {coup_object:.4e}")
    if separation_mas is not None:
        print(f"Separation: {separation_mas:.0f} mas")
    print(f"Raw Spectrum CSV:  {csv_outfile}")
    print(f"Output CSV:        {fiber_outfile}")
    print(f"Output Plot:       {png_outfile}")
    print("============================\n")

    fiber_flux = CubicSpline(wavelength_grid_angstroms, spec_star_ph)(wavelength_grid_angstroms) * coup_object

    with open(fiber_outfile, "w") as f:
        for j in range(len(wave)):
            f.write(f"{wave[j]:.3f},{fiber_flux[j]:.6e}\n")

    plt.figure(figsize=(8, 5))
    plt.plot(wave, fiber_flux[:len(wave)], label=f"Coupling={coup_str}")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Flux (photons/s)")
    plt.ticklabel_format(style='plain')
    #plt.yscale("log")
    plt.title(f"Spectrum for {filename_base} ({coup_str})")
    plt.legend()
    plt.tight_layout()

    # Restrict x-axis for non-star objects
    if "star" not in filename.lower():
        plt.xlim(655.5, 657.0)

    plt.savefig(png_outfile, dpi=150)
    plt.close()

    print(f"Saved -> {fiber_outfile} and {png_outfile}")

print("\nAll done!\n")
