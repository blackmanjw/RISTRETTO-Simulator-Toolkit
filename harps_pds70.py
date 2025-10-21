import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import os

# ---------------------------------
# File paths
# ---------------------------------
input_dir = 'input'
csv_file = os.path.join(input_dir, 'pds70_harps.csv')

txt_files = [
    'PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000_orig.txt',
    '2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000_orig.txt',
    'WISPIT2_star_BT-Settl-CIFIST-4400K-4logg_140000_orig.txt'
]

# ---------------------------------
# Safety checks
# ---------------------------------
if not os.path.exists(csv_file):
    raise FileNotFoundError(f"HARPS CSV file not found: {csv_file}")

# Load HARPS data once
data_csv = np.loadtxt(csv_file, delimiter=',', skiprows=1)
wavelength_csv = data_csv[:, 0]
flux_csv = data_csv[:, 1]

# ---------------------------------
# Colors for plotting
# ---------------------------------
colors_model = ['red', 'blue', 'orange']
colors_combined = ['purple', 'green', 'brown']

# ---------------------------------
# Loop through all TXT files
# ---------------------------------
for i, txt_name in enumerate(txt_files):
    txt_file = os.path.join(input_dir, txt_name)
    
    if not os.path.exists(txt_file):
        print(f"⚠️ Model TXT file not found, skipping: {txt_file}")
        continue

    data_txt = np.loadtxt(txt_file, comments='#')
    wavelength_txt = data_txt[:, 0]
    flux_txt = data_txt[:, 1]

    # ---------------------------------
    # Determine scaling factor near Hα
    # ---------------------------------
    region_mask_model = (wavelength_txt > 6550) & (wavelength_txt < 6557)
    region_mask_harps = (wavelength_csv > 6550) & (wavelength_csv < 6557)
    if np.any(region_mask_model) and np.any(region_mask_harps):
        scale_factor = np.median(flux_txt[region_mask_model]) / np.median(flux_csv[region_mask_harps])
    else:
        scale_factor = 0.025  # fallback
    print(f"{txt_name}: Using HARPS scaling factor: {scale_factor:.4f}")

    # ---------------------------------
    # Replace Hα region (6557–6565.6 Å) in BT-Settl with HARPS data
    # ---------------------------------
    transition_start = 6557
    transition_end = 6565.6

    harps_region_mask = (wavelength_csv >= transition_start) & (wavelength_csv <= transition_end)
    wavelength_harps_region = wavelength_csv[harps_region_mask]
    flux_harps_region = flux_csv[harps_region_mask] * scale_factor

    harps_interp = np.interp(
        wavelength_txt,
        wavelength_harps_region,
        flux_harps_region,
        left=np.nan,
        right=np.nan
    )

    flux_combined = flux_txt.copy()

    replace_start = 6558.0
    replace_end = 6565.4

    blend_mask = (wavelength_txt >= transition_start) & (wavelength_txt < replace_start)
    replace_mask = (wavelength_txt >= replace_start) & (wavelength_txt <= replace_end)
    blend_out_mask = (wavelength_txt > replace_end) & (wavelength_txt <= transition_end)

    if np.any(blend_mask):
        weights_in = (wavelength_txt[blend_mask] - transition_start) / (replace_start - transition_start)
        flux_combined[blend_mask] = (1 - weights_in) * flux_txt[blend_mask] + weights_in * harps_interp[blend_mask]
    if np.any(replace_mask):
        flux_combined[replace_mask] = harps_interp[replace_mask]
    if np.any(blend_out_mask):
        weights_out = (wavelength_txt[blend_out_mask] - replace_end) / (transition_end - replace_end)
        flux_combined[blend_out_mask] = (1 - weights_out) * harps_interp[blend_out_mask] + weights_out * flux_txt[blend_out_mask]

    # ---------------------------------
    # Define output TXT filename (drop '_orig')
    # ---------------------------------
    base_name = os.path.basename(txt_file)
    plot_base_name = os.path.splitext(base_name)[0].replace('_orig', '')
    output_txt = os.path.join(input_dir, f"{plot_base_name}.txt")

    # Custom header
    header_text = (
        f"# {plot_base_name} plus H-alpha peak from HARPS spectrum ({os.path.basename(csv_file)})\n"
        f"# Wavelength  Flux\n"
        f"# HARPS scaled by {scale_factor:.4f}\n"
        f"# Blended region: {replace_start}-{replace_end} Å"
    )

    # Save combined TXT
    np.savetxt(
        output_txt,
        np.column_stack([wavelength_txt, flux_combined]),
        fmt='%.6e',
        header=header_text,
        comments=''
    )

    print(f"✅ Output TXT created: {output_txt}")

    # ---------------------------------
    # Plot 1: Zoomed Hα region (6500–6650 Å)
    # ---------------------------------
    x_min, x_max = 6500, 6650
    plt.figure(figsize=(10, 6))
    plt.plot(wavelength_csv, flux_csv, label='HARPS Data', color='blue', alpha=0.5)
    plt.plot(wavelength_csv, flux_csv * scale_factor, label=f'HARPS × {scale_factor:.3f}', color='green', alpha=0.5)
    plt.plot(wavelength_txt, flux_txt, label='BT-Settl Model', color='red', alpha=0.4)
    plt.plot(wavelength_txt, flux_combined, label='Combined Spectrum', color='purple', linewidth=2)
    plt.axvspan(transition_start, transition_end, color='grey', alpha=0.2, label='Blended Region')
    plt.xlabel('Wavelength [Å]')
    plt.ylabel('Flux')
    plt.xlim(x_min, x_max)
    plt.yscale('log')
    
    # Auto-scale y-limits based on data in range
    mask_combined = (wavelength_txt >= x_min) & (wavelength_txt <= x_max)
    mask_harps = (wavelength_csv >= x_min) & (wavelength_csv <= x_max)
    all_flux = np.concatenate([flux_combined[mask_combined], flux_txt[mask_combined], flux_csv[mask_harps]*scale_factor])
    plt.ylim(all_flux.min()*0.9, all_flux.max()*1.1)
    
    plt.title(f'{plot_base_name} - Hα Zoom')
    plt.legend()
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(input_dir, f"{plot_base_name}_zoomed.png"), dpi=300)
    plt.close()
    
    # ---------------------------------
    # Plot 2: Extended range (6200–8600 Å)
    # ---------------------------------
    x_min, x_max = 6200, 8600
    plt.figure(figsize=(12, 6))
    plt.plot(wavelength_txt, flux_txt, color='red', alpha=0.3, label='BT-Settl Model')
    plt.plot(wavelength_txt, flux_combined, color='purple', linewidth=1.5, label='Combined Spectrum')
    plt.axvspan(transition_start, transition_end, color='grey', alpha=0.2, label='Hα Blended Region')
    plt.xlabel('Wavelength [Å]')
    plt.ylabel('Flux')
    plt.xlim(x_min, x_max)
    plt.yscale('log')
    
    mask_combined = (wavelength_txt >= x_min) & (wavelength_txt <= x_max)
    all_flux = np.concatenate([flux_combined[mask_combined], flux_txt[mask_combined]])
    plt.ylim(all_flux.min()*0.9, all_flux.max()*1.1)
    
    plt.title(f'{plot_base_name} - Full Spectrum')
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(input_dir, f"{plot_base_name}_full.png"), dpi=300)
    plt.close()
    
    # ---------------------------------
    # Plot 3: Hα Super-Zoom (6550–6575 Å)
    # ---------------------------------
    x_min, x_max = 6550, 6575
    plt.figure(figsize=(12, 6))
    plt.plot(wavelength_txt, flux_txt, color='red', alpha=0.3, label='BT-Settl Model')
    plt.plot(wavelength_txt, flux_combined, color='purple', linewidth=1.5, label='Combined Spectrum')
    plt.plot(wavelength_csv, flux_csv * scale_factor, color='green', alpha=0.5, linestyle='--', label=f'HARPS × {scale_factor:.3f}')
    plt.axvspan(transition_start, transition_end, color='grey', alpha=0.2, label='Hα Blended Region')
    plt.xlabel('Wavelength [Å]')
    plt.ylabel('Flux')
    plt.xlim(x_min, x_max)
    plt.yscale('log')
    
    mask_combined = (wavelength_txt >= x_min) & (wavelength_txt <= x_max)
    mask_harps = (wavelength_csv >= x_min) & (wavelength_csv <= x_max)
    all_flux = np.concatenate([flux_combined[mask_combined], flux_txt[mask_combined], flux_csv[mask_harps]*scale_factor])
    plt.ylim(all_flux.min()*0.9, all_flux.max()*1.1)
    
    plt.title(f'{plot_base_name} - Hα Super-Zoom')
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(input_dir, f"{plot_base_name}_superzoom.png"), dpi=300)
    plt.close()


    print(f"✅ Plots created for: {plot_base_name}\n")

# ---------------------------------
# Comparison figure across all models (Hα zoom)
# ---------------------------------
x_min, x_max = 6500, 6575
plt.figure(figsize=(14, 6))

all_fluxes = []  # will collect all flux values in the plotted range

for i, txt_name in enumerate(txt_files):
    txt_file = os.path.join(input_dir, txt_name)
    if not os.path.exists(txt_file):
        continue

    data_txt = np.loadtxt(txt_file, comments='#')
    wavelength_txt = data_txt[:, 0]
    flux_txt = data_txt[:, 1]

    # Determine scaling factor near Hα
    region_mask_model = (wavelength_txt > 6550) & (wavelength_txt < 6557)
    region_mask_harps = (wavelength_csv > 6550) & (wavelength_csv < 6557)
    if np.any(region_mask_model) and np.any(region_mask_harps):
        scale_factor = np.median(flux_txt[region_mask_model]) / np.median(flux_csv[region_mask_harps])
    else:
        scale_factor = 0.025

    # Replace Hα region
    transition_start = 6557
    transition_end = 6565.6
    harps_region_mask = (wavelength_csv >= transition_start) & (wavelength_csv <= transition_end)
    wavelength_harps_region = wavelength_csv[harps_region_mask]
    flux_harps_region = flux_csv[harps_region_mask] * scale_factor
    harps_interp = np.interp(wavelength_txt, wavelength_harps_region, flux_harps_region, left=np.nan, right=np.nan)

    flux_combined = flux_txt.copy()
    replace_start = 6558.0
    replace_end = 6565.4
    blend_mask = (wavelength_txt >= transition_start) & (wavelength_txt < replace_start)
    replace_mask = (wavelength_txt >= replace_start) & (wavelength_txt <= replace_end)
    blend_out_mask = (wavelength_txt > replace_end) & (wavelength_txt <= transition_end)

    if np.any(blend_mask):
        weights_in = (wavelength_txt[blend_mask] - transition_start) / (replace_start - transition_start)
        flux_combined[blend_mask] = (1 - weights_in) * flux_txt[blend_mask] + weights_in * harps_interp[blend_mask]
    if np.any(replace_mask):
        flux_combined[replace_mask] = harps_interp[replace_mask]
    if np.any(blend_out_mask):
        weights_out = (wavelength_txt[blend_out_mask] - replace_end) / (transition_end - replace_end)
        flux_combined[blend_out_mask] = (1 - weights_out) * harps_interp[blend_out_mask] + weights_out * flux_txt[blend_out_mask]

    # Clean base name
    base_name = os.path.basename(txt_file)
    plot_base_name = os.path.splitext(base_name)[0].replace('_orig', '')

    # Plot original model and combined
    plt.plot(wavelength_txt, flux_txt, color=colors_model[i], alpha=0.4, label=f'{plot_base_name} (model)')
    plt.plot(wavelength_txt, flux_combined, color=colors_combined[i], linewidth=1.5, label=f'{plot_base_name} (Hα blended)')

    # Collect flux values in the x-range for auto-scaling
    mask_combined = (wavelength_txt >= x_min) & (wavelength_txt <= x_max)
    mask_harps = (wavelength_csv >= x_min) & (wavelength_csv <= x_max)
    all_fluxes.append(flux_txt[mask_combined])
    all_fluxes.append(flux_combined[mask_combined])
    all_fluxes.append(flux_csv[mask_harps] * scale_factor)

# Concatenate all flux values in range
all_fluxes = np.concatenate(all_fluxes)

# Plot scaled HARPS
plt.plot(wavelength_csv, flux_csv * scale_factor, color='black', linestyle='--', label='HARPS scaled', alpha=0.5)

plt.xlabel('Wavelength [Å]')
plt.ylabel('Flux')
plt.title('Comparison of Models with Hα Blending')
plt.xlim(x_min, x_max)
plt.yscale('log')

# Auto-scale y-axis with 10% padding
plt.ylim(all_fluxes.min() * 0.9, all_fluxes.max() * 1.1)

plt.grid(True, which='both', ls='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(input_dir, "comparison_Ha_blended.png"), dpi=300)
plt.show()

print("✅ Comparison figure created: comparison_Ha_blended.png")
