#!/usr/bin/env python3
from pyechelle.simulator import Simulator
from pyechelle.spectrograph import ZEMAX
from pyechelle.sources import CSVSource
import numpy as np
import matplotlib.pyplot as plt
import smplotlib  # ✅ for consistent style
import os
from matplotlib.ticker import ScalarFormatter  # ✅ for plain number formatting

# === Common parameters ===
base_path = '/storage/homefs/jb23l046/Simu_run/data/'
exp_time = 3600  # seconds
spectrograph_model = "RISTRETTOtest18"

# === Define systems and their planets ===
systems = {
    "2MJ1612": {
        "star_base": "2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000",
        "planets": ["2MJ1612b_BT-Settl-CIFIST-1200K-3.5logg_140000"],
    },
    "WISPIT2": {
        "star_base": "WISPIT2_star_BT-Settl-CIFIST-4400K-4logg_140000",
        "planets": ["WISPIT2b_BT-Settl-CIFIST-1400K-4logg_140000"],
    },
    "PDS70": {
        "star_base": "PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000",
        "planets": [
            "PDS70b_BT-Settl-CIFIST-1400K-4logg_140000",
            "PDS70c_BT-Settl-CIFIST-1300K-4logg_140000",
        ],
    },
}

# === Variations for Ha lines ===
ha_variations = [
    ("60", "12"),
    ("60", "14"),
    ("80", "12"),
    ("80", "14"),
]

# === Helper functions ===
def planet_filename(base, ha1, ha2):
    return f"{base}_Ha_{ha1}_{ha2}_0.45.csv"

def extract_ha_tag(filename):
    """Extract Ha tag (e.g., Ha_80_12) without trailing numeric parts."""
    parts = filename.replace(".csv", "").split("_")
    try:
        ha_idx = parts.index("Ha")
        return "_".join(parts[ha_idx:ha_idx+3])
    except ValueError:
        return "Ha_unknown"

def disable_sci_notation(ax):
    """Disable scientific notation on x-axis."""
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.xaxis.get_major_formatter().set_useOffset(False)

# === Main Loop ===
for system_name, config in systems.items():
    print(f"\n=== Processing system: {system_name} ===")

    star_base = config["star_base"]
    star_file_045 = f"{star_base}_0.45.csv"
    star_file_350e4 = f"{star_base}_3.50e-04.csv"

    # Load the secondary star file once per system
    star_350e4_data = np.loadtxt(os.path.join(base_path, "in", star_file_350e4), delimiter=",")

    for planet_base in config["planets"]:
        for ha1, ha2 in ha_variations:
            planet_file = planet_filename(planet_base, ha1, ha2)
            planet_path = os.path.join(base_path, "in", planet_file)

            # --- Load planet data ---
            planet_data = np.loadtxt(planet_path, delimiter=",")

            # --- Check wavelength match ---
            if not np.allclose(star_350e4_data[:, 0], planet_data[:, 0]):
                raise ValueError(f"Wavelength grids mismatch between {planet_file} and {star_file_350e4}")

            # --- Combine star + planet flux ---
            combined_flux = planet_data[:, 1] + star_350e4_data[:, 1]
            combined_data = np.column_stack((planet_data[:, 0], combined_flux))

            # --- Build output filenames ---
            ha_tag_str = extract_ha_tag(planet_file)
            planet_name = os.path.basename(planet_base.split('_')[0])  # e.g. "WISPIT2b"
            combined_basename = f"{planet_name}_star_plus_planet_{ha_tag_str}"
            combined_file = f"{combined_basename}.csv"
            combined_path = os.path.join(base_path, "in", combined_file)

            np.savetxt(combined_path, combined_data, delimiter=",", fmt="%.6e")
            print(f"  → Combined CSV saved: {combined_path}")

            # === ✅ Plot and save with same filename ===
            wavelength = combined_data[:, 0]
            flux = combined_data[:, 1]

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(wavelength, flux, color="black", linewidth=1.0)
            ax.set_xlabel("Wavelength (Å)")
            ax.set_ylabel(r"Flux (ph s$^{-1}$ Å$^{-1}$)")
            ax.set_title(f"{planet_name} {ha_tag_str} - Star + Planet Spectrum")
            ax.set_yscale("log")
            disable_sci_notation(ax)
            fig.tight_layout()

            combined_plot = os.path.join(base_path, "in", f"{combined_basename}.png")
            fig.savefig(combined_plot, dpi=300)
            plt.close(fig)
            print(f"  📊 Plot saved: {combined_plot}")

            # === ✅ Create a second plot showing both star-only and combined spectra ===
            star_wavelength = star_350e4_data[:, 0]
            star_flux = star_350e4_data[:, 1]

            fig2, ax2 = plt.subplots(figsize=(10, 6))
            ax2.plot(star_wavelength, star_flux, color="tab:red", linewidth=1.0, label="Star only")
            ax2.plot(wavelength, flux, color="black", linewidth=1.0, label="Star + Planet")

            ax2.set_xlabel("Wavelength (Å)")
            ax2.set_ylabel(r"Flux (ph s$^{-1}$ Å$^{-1}$)")
            ax2.set_title(f"{planet_name} {ha_tag_str} - Star vs Star+Planet Spectrum")
            ax2.set_yscale("log")
            disable_sci_notation(ax2)
            ax2.legend()
            fig2.tight_layout()

            combined_plot_with_star = os.path.join(base_path, "in", f"{combined_basename}_combined.png")
            fig2.savefig(combined_plot_with_star, dpi=300)
            plt.close(fig2)
            print(f"  📊 Combined plot saved: {combined_plot_with_star}")


            # --- Run simulator ---
            sim = Simulator(ZEMAX(spectrograph_model))
            sim.set_ccd(1)
            sim.set_fibers([1, 2])
            sim.set_sources([
                CSVSource(os.path.join(base_path, "in", star_file_045), list_like=False,
                          wavelength_units='nm', flux_units='ph/s/AA'),
                CSVSource(combined_path, list_like=False,
                          wavelength_units='nm', flux_units='ph/s/AA')
            ])
            sim.set_exposure_time(exp_time)
            sim.set_bias(250)
            sim.set_read_noise(3)

            # --- Output FITS filename ---
            output_filename = f"{combined_basename}_{exp_time}s.fits"
            output_path = os.path.join(base_path, "out/onaxis", output_filename)

            sim.set_output(output_path, overwrite=True)
            sim.run()
            print(f"  ✓ Simulation complete for {planet_file} → {output_filename}")

print("\n=== All simulations complete! ===")
