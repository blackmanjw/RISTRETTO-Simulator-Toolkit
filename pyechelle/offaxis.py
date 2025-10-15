from pyechelle.simulator import Simulator
from pyechelle.spectrograph import ZEMAX
from pyechelle.sources import CSVSource
import numpy as np
import os
import re

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

# === Planet Ha variations ===
ha_variations = [
    ("60", "12"),
    ("60", "14"),
    ("80", "12"),
    ("80", "14"),
]

# === Star-only 0.45 CSV files ===
stars_045 = {
    "2MJ1612": "2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000_0.45.csv",
    "WISPIT2": "WISPIT2_star_BT-Settl-CIFIST-4400K-4logg_140000_0.45.csv",
    "PDS70":  "PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000_0.45.csv",
}

# === Directories ===
in_dir = os.path.join(base_path, "in")
in_offaxis_dir = os.path.join(in_dir, "offaxis")
out_dir = os.path.join(base_path, "out/offaxis")

# Ensure subdirectories exist
os.makedirs(in_offaxis_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

# === Helper functions ===
def extract_mas_tag(filename):
    match = re.search(r'_(\d{3})mas', filename)
    return match.group(1) + "mas" if match else "unknown_mas"

def extract_coupling(filename):
    """Extract the coupling value like 1.60e-04 from filename."""
    match = re.search(r'_([0-9]\.\d+e-\d{2})_', filename)
    return match.group(1) if match else "unknown_coupling"

def extract_ha_tag(filename):
    parts = filename.replace(".csv", "").split("_")
    try:
        ha_idx = parts.index("Ha")
        return "_".join(parts[ha_idx:ha_idx+3])
    except ValueError:
        return "Ha_unknown"

# ==============================
# Part 1: Star-only 0.45 simulations
# ==============================
print("\n=== Running star-only 0.45 simulations ===")
for star_name, star_file in stars_045.items():
    star_path = os.path.join(in_dir, star_file)
    if not os.path.exists(star_path):
        print(f"⚠️ Star file not found: {star_file}")
        continue

    print(f"\nRunning star-only simulation: {star_file}")
    sim = Simulator(ZEMAX(spectrograph_model))
    sim.set_ccd(1)
    sim.set_fibers([1])
    sim.set_sources([
        CSVSource(star_path, list_like=False,
                  wavelength_units='nm', flux_units='ph/s/AA')
    ])
    sim.set_exposure_time(exp_time)
    sim.set_bias(250)
    sim.set_read_noise(3)

    output_filename = f"{star_name}_star_0.45_{exp_time}s.fits"
    output_path = os.path.join(out_dir, output_filename)

    sim.set_output(output_path, overwrite=True)
    sim.run()
    print(f"✓ Star-only simulation complete → {output_filename}")

# ==============================
# Part 2: Off-axis star + planet simulations
# ==============================
print("\n=== Running off-axis star+planet simulations ===")

for system_name, config in systems.items():
    print(f"\n==============================")
    print(f"=== Processing system: {system_name} ===")
    print(f"==============================")

    star_base = config["star_base"]
    star_pattern = f"{star_base}_"

    # Find all off-axis star CSVs
    star_files = [
        f for f in os.listdir(in_dir)
        if f.startswith(star_pattern) and re.search(r'_(\d\.\d+e-\d{2})_\d{3}mas\.csv', f)
    ]
    if not star_files:
        print(f"⚠️ No off-axis star files found for {system_name}")
        continue

    for star_file in sorted(star_files):
        star_path = os.path.join(in_dir, star_file)
        mas_tag = extract_mas_tag(star_file)
        coupling = extract_coupling(star_file)
        print(f"\n--- Using star {star_file} ({coupling}, {mas_tag}) ---")

        star_data = np.loadtxt(star_path, delimiter=",")

        for planet_base in config["planets"]:
            planet_name = os.path.basename(planet_base.split('_')[0])

            for ha1, ha2 in ha_variations:
                planet_file = f"{planet_base}_Ha_{ha1}_{ha2}_0.45.csv"
                planet_path = os.path.join(in_dir, planet_file)

                if not os.path.exists(planet_path):
                    print(f"⚠️ Missing planet file: {planet_file}")
                    continue

                planet_data = np.loadtxt(planet_path, delimiter=",")

                # Check wavelength alignment
                if not np.allclose(star_data[:, 0], planet_data[:, 0]):
                    raise ValueError(f"Wavelength mismatch: {planet_file} and {star_file}")

                # Combine fluxes
                combined_flux = planet_data[:, 1] + star_data[:, 1]
                combined_data = np.column_stack((planet_data[:, 0], combined_flux))

                ha_tag = extract_ha_tag(planet_file)
                combined_file = f"{system_name}_star_plus_planet_{ha_tag}_{coupling}_{mas_tag}.csv"
                combined_path = os.path.join(in_offaxis_dir, combined_file)
                np.savetxt(combined_path, combined_data, delimiter=",", fmt="%.6e")
                print(f"  → Combined spectrum saved: {combined_file}")

                # --- Run simulator ---
                sim = Simulator(ZEMAX(spectrograph_model))
                sim.set_ccd(1)
                sim.set_fibers([1])
                sim.set_sources([
                    CSVSource(combined_path, list_like=False,
                              wavelength_units='nm', flux_units='ph/s/AA')
                ])
                sim.set_exposure_time(exp_time)
                sim.set_bias(250)
                sim.set_read_noise(3)

                # --- Output FITS file ---
                output_filename = f"{planet_name}_{ha_tag}_{coupling}_{mas_tag}_{exp_time}s.fits"
                output_path = os.path.join(out_dir, output_filename)

                sim.set_output(output_path, overwrite=True)
                sim.run()
                print(f"  ✓ Simulation complete → {output_filename}")

print("\n=== ✅ All simulations complete! ===")
