import numpy as np
import os
import matplotlib.pyplot as plt

# =========================
# User settings
# =========================
plot_flag = True          # Set True to plot each file
use_trapz = False         # Set True to use np.trapz, False to use sum()/dl method

# Constants
pc_to_cm = 3.086e18
R_jupiter = 7.1492e9  # cm

# Files and stars
input_files = ['Ha_60_14.dat', 'Ha_60_12.dat', 'Ha_80_14.dat', 'Ha_80_12.dat']

stars = {
    'pds70b': {'d_pc': 113.4, 'R_star': 2.0 * R_jupiter,'halpha_flux': 8.1e-16},
    'pds70c': {'d_pc': 113.4, 'R_star': 1.6 * R_jupiter,'halpha_flux': 3.1e-16},
    'WISPIT2': {'d_pc': 133.0, 'R_star': 1.6 * R_jupiter,'halpha_flux': 1.38e-15},
    '2MJ1612b': {'d_pc': 131.9, 'R_star': 1.5 * R_jupiter,'halpha_flux': 8.2e-16}
}

# =========================
# Processing loop
# =========================
for star_name, params in stars.items():
    d_cm = params['d_pc'] * pc_to_cm
    R_star = params['R_star']
    target_line_flux = params['halpha_flux'] # erg/s/cm2
    
    # Create output directory
    out_dir = f'Halpha_model/scaled/{star_name}'
    os.makedirs(out_dir, exist_ok=True)
    
    for filename in input_files:
        # Load data
        filepath = f'Halpha_model/{filename}'
        data = np.loadtxt(filepath, comments='#')
        wavelength = data[:,0]
        Fsurf = data[:,1]
        
        # Step 0: ensure wavelengths increasing for scaling
        if wavelength[0] > wavelength[-1]:
            wavelength = wavelength[::-1]
            Fsurf = Fsurf[::-1]
        
        # Step 1: scale to Earth distance and multiply by 0.01
        F_earth = Fsurf * (R_star/d_cm)**2 * 0.01
        
        # Step 2: scale to target line flux
        if use_trapz:
            # Trapezoidal integration method
            F_line = np.trapz(F_earth, wavelength)
            scale_factor = target_line_flux / F_line
            F_earth_scaled = F_earth * scale_factor
        else:
            # sum()/dl method
            dl = wavelength[1] - wavelength[0]  # assumes uniform spacing
            F_earth_scaled = F_earth / F_earth.sum() / dl * target_line_flux
        
        # Step 3: convert back to surface flux
        Fsurf_new = F_earth_scaled * (d_cm / R_star)**2
        
        # Step 4: reverse to original decreasing wavelength if needed
        if data[0,0] > data[-1,0]:
            wavelength_out = wavelength[::-1]
            Fsurf_new_out = Fsurf_new[::-1]
        else:
            wavelength_out = wavelength
            Fsurf_new_out = Fsurf_new
        
        # Step 5: save to .dat with original formatting
        out_path = os.path.join(out_dir, filename)
        with open(out_path, 'w') as f:
            f.write('#Wavelength[A]\tEnergy Flux Density[erg/s/A/cm2]\n')
            for wl, flux in zip(wavelength_out, Fsurf_new_out):
                f.write(f'{wl:20.6f}\t{flux:20.10e}\n')
        
        # Print info
        method_name = "trapz" if use_trapz else "sum()/dl"
        print(f'{star_name}/{filename} processed using {method_name}.')
        
        # Optional plotting
        if plot_flag:
            plt.figure(figsize=(10,6))
            plt.plot(wavelength_out, Fsurf[::-1] if data[0,0] > data[-1,0] else Fsurf, 
                     label='Original Surface Flux')
            plt.plot(wavelength_out, Fsurf_new_out, linestyle='--', label='Scaled Surface Flux')
            plt.xlabel('Wavelength [Å]')
            plt.ylabel('Flux [erg/s/Å/cm²]')
            plt.title(f'{star_name}: {filename} ({method_name})')
            plt.yscale('log')
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()
