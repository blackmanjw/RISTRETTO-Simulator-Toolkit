import numpy as np
import matplotlib.pyplot as plt
import os
import smplotlib


# --- Constants ---
R_sun = 6.957e10        # cm
R_jup = 7.1492e9        # cm
pc_to_cm = 3.0857e18    # cm

# --- Input Parameters ---
distance_pc = 113
R_star = 1.26 * R_sun
R_planet = 1.6 * R_jup
distance_cm = distance_pc * pc_to_cm

# --- Ensure output directory exists ---
os.makedirs("pds70", exist_ok=True)

# --- Read Data (surface spectra) ---
star_data = np.loadtxt("input/PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt", comments="#")
planet_data = np.loadtxt("input/PDS70b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12.txt", comments="#")

# --- Extract columns ---
wl_star, flux_star_surface = star_data[:, 0], star_data[:, 1]
wl_planet, flux_planet_surface = planet_data[:, 0], planet_data[:, 1]

# --- Scale flux to Earth's distance ---
flux_star_earth = flux_star_surface * (R_star / distance_cm)**2
flux_planet_earth = flux_planet_surface * (R_planet / distance_cm)**2

# --- Plot range ---
x_min, x_max = 6200, 8500

def restrict_range(wl, flux, x_min, x_max):
    mask = (wl >= x_min) & (wl <= x_max)
    return wl[mask], flux[mask]

# Restrict data
wl_star_sub, flux_star_surface_sub = restrict_range(wl_star, flux_star_surface, x_min, x_max)
wl_planet_sub, flux_planet_surface_sub = restrict_range(wl_planet, flux_planet_surface, x_min, x_max)
wl_star_earth_sub, flux_star_earth_sub = restrict_range(wl_star, flux_star_earth, x_min, x_max)
wl_planet_earth_sub, flux_planet_earth_sub = restrict_range(wl_planet, flux_planet_earth, x_min, x_max)

# --- Plot 1: Flux at the surface ---
plt.figure(figsize=(10, 6))
plt.plot(wl_star_sub, flux_star_surface_sub, label="PDS 70 (star surface)", color="goldenrod")
plt.plot(wl_planet_sub, flux_planet_surface_sub, label="PDS 70 b (planet surface)", color="mediumblue")
plt.xlim(x_min, x_max)
plt.yscale("log")
plt.xlabel("Wavelength (Å)", fontsize=12)
plt.ylabel(r"Flux ($\mathrm{erg\,cm^{-2}\,s^{-1}\,\AA^{-1}}$)", fontsize=12)
plt.title("PDS 70 and PDS 70b at Planet/Star Surface (BT-Settl) + Aoyama H-alpha (Ha_60_12)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.ticklabel_format(axis="x", style="plain")
plt.tight_layout()
plt.savefig("pds70/plot1_surface_flux.png", dpi=300)
plt.close()

# --- Plot 2: Flux at Earth ---
plt.figure(figsize=(10, 6))
plt.plot(wl_star_earth_sub, flux_star_earth_sub, label="PDS 70 (at Earth, 113 pc)", color="goldenrod")
plt.plot(wl_planet_earth_sub, flux_planet_earth_sub, label="PDS 70 b (at Earth, 113 pc)", color="mediumblue")
plt.xlim(x_min, x_max)
plt.yscale("log")
plt.xlabel("Wavelength (Å)", fontsize=12)
plt.ylabel(r"Flux ($\mathrm{erg\,cm^{-2}\,s^{-1}\,\AA^{-1}}$)", fontsize=12)
plt.title("PDS 70 and PDS 70b at Earth (BT-Settl) + Aoyama H-alpha (Ha_60_12)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.ticklabel_format(axis="x", style="plain")
plt.tight_layout()
plt.savefig("pds70/plot2_earth_flux.png", dpi=300)
plt.close()

# --- Plot 3: Direct CSV flux comparison (already in nm) ---
files = {
    "PDS70b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12_0.45.csv": "External Spaxel - Planet (0.45)",
    "PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000_3.50e-04.csv": "External Spaxel - Star (3.50e-04)",
    "PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000_0.45.csv": "Centre Spaxel - Star (0.45)",
}

plt.figure(figsize=(10, 6))

for filename, label in files.items():
    path = os.path.join("output", filename)
    data = np.loadtxt(path, delimiter=",", comments="#")
    wl_nm = data[:, 0]  # already in nm
    flux = data[:, 1]   # use as-is
    plt.plot(wl_nm, flux, label=label)

plt.xscale("linear")
plt.yscale("log")
plt.xlabel("Wavelength (nm)", fontsize=12)
plt.ylabel("Flux (photons/s) at telescope", fontsize=12)
plt.title("PDS 70b observation - Pyechelle Input")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("pds70/plot3_csv_flux_comparison.png", dpi=300)
plt.close()

print("✅ All three plots saved in the 'pds70/' directory.")
