import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import os

# ---------------------------------
# File paths
# ---------------------------------
input_dir = 'input'
csv_file = os.path.join(input_dir, 'pds70_harps.csv')
txt_file = os.path.join(input_dir, 'PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt')

# ---------------------------------
# Safety checks
# ---------------------------------
if not os.path.exists(csv_file):
    raise FileNotFoundError(f"HARPS CSV file not found: {csv_file}")
if not os.path.exists(txt_file):
    raise FileNotFoundError(f"Model TXT file not found: {txt_file}")

# ---------------------------------
# Load data
# ---------------------------------
data_csv = np.loadtxt(csv_file, delimiter=',', skiprows=1)
wavelength_csv = data_csv[:, 0]
flux_csv = data_csv[:, 1]

data_txt = np.loadtxt(txt_file, comments='#')
wavelength_txt = data_txt[:, 0]
flux_txt = data_txt[:, 1]

# ---------------------------------
# Determine scaling factor automatically near Hα region
# ---------------------------------
region_mask_model = (wavelength_txt > 6550) & (wavelength_txt < 6557)
region_mask_harps = (wavelength_csv > 6550) & (wavelength_csv < 6557)
if np.any(region_mask_model) and np.any(region_mask_harps):
    scale_factor = np.median(flux_txt[region_mask_model]) / np.median(flux_csv[region_mask_harps])
else:
    scale_factor = 0.025  # fallback
print(f"Using HARPS scaling factor: {scale_factor:.4f}")

# ---------------------------------
# Replace Hα region (6557–6565.6 Å) in BT-Settl with HARPS data
# ---------------------------------
transition_start = 6557
transition_end = 6565.6

# HARPS region selection and scaling
harps_region_mask = (wavelength_csv >= transition_start) & (wavelength_csv <= transition_end)
wavelength_harps_region = wavelength_csv[harps_region_mask]
flux_harps_region = flux_csv[harps_region_mask] * scale_factor

# Interpolate HARPS flux to model wavelength grid
harps_interp = np.interp(
    wavelength_txt,
    wavelength_harps_region,
    flux_harps_region,
    left=np.nan,
    right=np.nan
)

# Combined flux array (copy model)
flux_combined = flux_txt.copy()

# ---------------------------------
# Replace 6558–6565.4 completely with HARPS data
# ---------------------------------
replace_start = 6558.0
replace_end = 6565.4

# Define masks
blend_mask = (wavelength_txt >= transition_start) & (wavelength_txt < replace_start)
replace_mask = (wavelength_txt >= replace_start) & (wavelength_txt <= replace_end)
blend_out_mask = (wavelength_txt > replace_end) & (wavelength_txt <= transition_end)

# Copy original model
flux_combined = flux_txt.copy()

# (1) Blend-in zone (6557–6558 Å)
if np.any(blend_mask):
    weights_in = (wavelength_txt[blend_mask] - transition_start) / (replace_start - transition_start)
    flux_combined[blend_mask] = (
        (1 - weights_in) * flux_txt[blend_mask] + weights_in * harps_interp[blend_mask]
    )

# (2) Full replacement zone (6558–6565.4 Å)
if np.any(replace_mask):
    flux_combined[replace_mask] = harps_interp[replace_mask]

# (3) Blend-out zone (6565.4–6565.6 Å)
if np.any(blend_out_mask):
    weights_out = (wavelength_txt[blend_out_mask] - replace_end) / (transition_end - replace_end)
    flux_combined[blend_out_mask] = (
        (1 - weights_out) * harps_interp[blend_out_mask] + weights_out * flux_txt[blend_out_mask]
    )

# ---------------------------------
# Define output filenames based on input TXT file
# ---------------------------------
base_name = os.path.splitext(os.path.basename(txt_file))[0]  # e.g., PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000
output_root = f"{base_name}_harps"

# TXT output
output_txt = os.path.join(input_dir, f"{output_root}.txt")
np.savetxt(
    output_txt,
    np.column_stack([wavelength_txt, flux_combined]),
    fmt='%.6e',  # scientific notation for readability
    header=f'# Wavelength  Flux\n# HARPS scaled by {scale_factor:.4f}\n# Blended region: {replace_start}-{replace_end} Å',
    comments=''
)

# Plot 1: Zoomed region
plt.figure(figsize=(10, 6))
plt.plot(wavelength_csv, flux_csv, label='HARPS Data', color='blue', alpha=0.5)
plt.plot(wavelength_csv, flux_csv * scale_factor, label=f'HARPS × {scale_factor:.3f}', color='green', alpha=0.5)
plt.plot(wavelength_txt, flux_txt, label='BT-Settl Model', color='red', alpha=0.4)
plt.plot(wavelength_txt, flux_combined, label='Combined (smoothed)', color='purple', linewidth=2)
plt.axvspan(transition_start, transition_end, color='grey', alpha=0.2, label='Blended Region')
plt.xlabel('Wavelength [Å]')
plt.ylabel('Flux')
plt.xlim(6500, 6650)
plt.yscale('log')
plt.ylim(6e5, 3e8)
plt.title('PDS70 Spectrum vs Model (Zoomed on Hα)')
ax = plt.gca()
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.ticklabel_format(useOffset=False, axis='x')
plt.legend()
plt.grid(True, which='both', ls='--', alpha=0.5)
plt.tight_layout()
plt.savefig(os.path.join(input_dir, f"{output_root}_zoomed.png"), dpi=300)
plt.show()

# Plot 2: Extended range (6200–8600 Å)
plt.figure(figsize=(12, 6))
plt.plot(wavelength_txt, flux_txt, color='red', alpha=0.3, label='BT-Settl Model')
plt.plot(wavelength_txt, flux_combined, color='purple', linewidth=1.5, label='Combined Spectrum')
plt.axvspan(transition_start, transition_end, color='grey', alpha=0.2, label='Hα Blended Region')
plt.xlabel('Wavelength [Å]')
plt.ylabel('Flux')
plt.title('PDS70 Combined Spectrum (6200–8600 Å)')
plt.yscale('log')
plt.xlim(6200, 8600)
plt.ylim(7e5, 4e6)
plt.grid(True, which='both', ls='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(input_dir, f"{output_root}_full.png"), dpi=300)
plt.show()

# Plot 3: Hα Super-Zoom (6550–6575 Å)
plt.figure(figsize=(12, 6))
plt.plot(wavelength_txt, flux_txt, color='red', alpha=0.3, label='BT-Settl Model')
plt.plot(wavelength_txt, flux_combined, color='purple', linewidth=1.5, label='Combined Spectrum')
plt.plot(wavelength_csv, flux_csv * scale_factor, label=f'HARPS × {scale_factor:.3f}', color='green', alpha=0.5)
plt.axvspan(transition_start, transition_end, color='grey', alpha=0.2, label='Hα Blended Region')
plt.xlabel('Wavelength [Å]')
plt.ylabel('Flux')
plt.title('PDS70 Combined Spectrum (Hα Super-Zoom)')
plt.yscale('log')
plt.xlim(6550, 6575)
plt.ylim(9e5, 4e6)
plt.grid(True, which='both', ls='--', alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(input_dir, f"{output_root}_superzoom.png"), dpi=300)
plt.show()

print("✅ Output files created:")
print(f"   - {output_root}.txt")
print(f"   - {output_root}_zoomed.png")
print(f"   - {output_root}_full.png")
print(f"   - {output_root}_superzoom.png")
