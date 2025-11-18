import numpy as np
import os
import argparse
import matplotlib.pyplot as plt

# =========================
# User settings
# =========================
plot_flag = True          # Set True to plot each file
show_plots = False        # Set True to display plots interactively with plt.show()
use_trapz = False         # Set True to use np.trapz, False to use sum()/dl method

# Constants
pc_to_cm = 3.086e18
R_jupiter = 7.1492e9  # cm

# Files and stars
input_files = ['Ha_60_14.dat', 'Ha_60_12.dat', 'Ha_80_14.dat', 'Ha_80_12.dat']
default_stars = {
    'pds70b': {'d_pc': 113.4, 'R_planet': 2.0 * R_jupiter, 'halpha_flux': 8.1e-16},
    'pds70c': {'d_pc': 113.4, 'R_planet': 1.6 * R_jupiter, 'halpha_flux': 3.1e-16},
    'WISPIT2': {'d_pc': 133.0, 'R_planet': 1.6 * R_jupiter, 'halpha_flux': 1.38e-15},
    '2MJ1612b': {'d_pc': 131.9, 'R_planet': 1.5 * R_jupiter, 'halpha_flux': 8.2e-16},
}

# =========================
# Command-line arguments
# =========================
parser = argparse.ArgumentParser(description="Scale Halpha models using optional flux inputs.")

parser.add_argument("--pds70b", type=float, default=None, help="Halpha flux for PDS 70 b")
parser.add_argument("--pds70c", type=float, default=None, help="Halpha flux for PDS 70 c")
parser.add_argument("--WISPIT2", type=float, default=None, help="Halpha flux for WISPIT2")
parser.add_argument("--2MJ1612b", type=float, default=None, help="Halpha flux for 2MJ1612b")

args = parser.parse_args()

# Build final star dictionary using defaults unless overridden
stars = {}
for name, params in default_stars.items():
    user_flux = getattr(args, name.replace("-", "_"), None)
    stars[name] = params.copy()
    if user_flux is not None:
        stars[name]['halpha_flux'] = user_flux
        print(f"Using command-line Hα flux for {name}: {user_flux:e}")
    else:
        print(f"Using default Hα flux for {name}: {params['halpha_flux']:e}")

# =========================
# Processing loop
# =========================
for star_name, params in stars.items():
    d_cm = params['d_pc'] * pc_to_cm
    R_planet = params['R_planet']
    target_line_flux = params['halpha_flux']  # erg/s/cm2
    
    # Create output directory
    out_dir = f'Halpha_model/scaled/{star_name}'
    os.makedirs(out_dir, exist_ok=True)
    
    for filename in input_files:
        # Load data
        filepath = f'Halpha_model/{filename}'
        data = np.loadtxt(filepath, comments='#')
        wavelength = data[:, 0]
        Fsurf = data[:, 1]
        
        # Ensure increasing wavelength
        if wavelength[0] > wavelength[-1]:
            wavelength = wavelength[::-1]
            Fsurf = Fsurf[::-1]
        
        # Step 1: scale to Earth distance
        F_earth = Fsurf * (R_planet / d_cm) ** 2 * 0.01
        
        # Step 2: integrate and rescale
        if use_trapz:
            F_line = np.trapz(F_earth, wavelength)
            scale_factor = target_line_flux / F_line
            F_earth_scaled = F_earth * scale_factor
        else:
            dl = wavelength[1] - wavelength[0]
            F_earth_scaled = F_earth / F_earth.sum() / dl * target_line_flux
        
        # Step 3: convert back to surface flux
        Fsurf_new = F_earth_scaled * (d_cm / R_planet) ** 2
        
        # Step 4: reverse if needed
        if data[0, 0] > data[-1, 0]:
            wavelength_out = wavelength[::-1]
            Fsurf_new_out = Fsurf_new[::-1]
        else:
            wavelength_out = wavelength
            Fsurf_new_out = Fsurf_new
        
        # Step 5: save file
        out_path = os.path.join(out_dir, filename)
        with open(out_path, 'w') as f:
            f.write('#Wavelength[A]\tEnergy Flux Density[erg/s/A/cm2]\n')
            for wl, flux in zip(wavelength_out, Fsurf_new_out):
                f.write(f"{wl:20.6f}\t{flux:20.10e}\n")
        
        # Plotting
        method = "trapz" if use_trapz else "sum()/dl"
        print(f"{star_name}/{filename} processed using {method}.")
        
        if plot_flag:
            plt.figure(figsize=(10, 6))
            plt.plot(
                wavelength_out,
                Fsurf[::-1] if data[0, 0] > data[-1, 0] else Fsurf,
                label="Original Surface Flux"
            )
            plt.plot(wavelength_out, Fsurf_new_out, "--", label="Scaled Surface Flux")
            plt.xlabel("Wavelength [Å]")
            plt.ylabel("Flux [erg/s/Å/cm²]")
            plt.title(f"{star_name}: {filename} ({method})")
            plt.yscale("log")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            if show_plots:
                plt.show()
            else:
                plt.close()
