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
