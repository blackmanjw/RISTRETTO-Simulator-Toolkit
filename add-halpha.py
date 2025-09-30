# Run like this: python add-halpha.py 2MJ1612b_BT-Settl-CIFIST-1200K-3.5logg_140000.txt Ha_60_12.dat
#
#!/usr/bin/env python3
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) != 3:
        print("Usage: python combine_spectra.py <cifist_filename> <halpha_filename>")
        sys.exit(1)

    cifist_name = sys.argv[1]
    halpha_name = sys.argv[2]

    # Hard-coded directory structure
    cifist_path = os.path.join("input", cifist_name)
    halpha_path = os.path.join("Halpha_model", halpha_name)

    # Read spectra
    wl1, flux1 = np.loadtxt(cifist_path, unpack=True)
    wl2, flux2 = np.loadtxt(halpha_path, unpack=True)

    # Interpolate cifist flux onto halpha wavelength grid (for comparison)
    flux1_interp_on_halpha = np.interp(wl2, wl1, flux1)

    # Keep only H-alpha points where H-alpha flux > cifist flux
    mask_keep_halpha = flux2 > flux1_interp_on_halpha
    wl2_kept = wl2[mask_keep_halpha]
    flux2_kept = flux2[mask_keep_halpha]

    # Remove cifist points in H-alpha wavelength range
    wl2_min = wl2_kept.min() if len(wl2_kept) > 0 else wl2.min()
    wl2_max = wl2_kept.max() if len(wl2_kept) > 0 else wl2.max()
    mask_keep_cifist = (wl1 < wl2_min) | (wl1 > wl2_max)
    wl1_clean = wl1[mask_keep_cifist]
    flux1_clean = flux1[mask_keep_cifist]

    # Combine the spectra
    wl_combined = np.concatenate([wl1_clean, wl2_kept])
    flux_combined = np.concatenate([flux1_clean, flux2_kept])

    # Sort by wavelength
    sorted_idx = np.argsort(wl_combined)
    wl_combined = wl_combined[sorted_idx]
    flux_combined = flux_combined[sorted_idx]

    # Output filenames
    base1 = os.path.splitext(cifist_name)[0]
    base2 = os.path.splitext(halpha_name)[0]
    out_file = f"{base1}_{base2}.dat"
    combined_plot = f"{base1}_{base2}.png"
    halpha_plot = f"{base2}_input.png"
    inputs_plot = f"{base1}_{base2}_inputs.png"

    # Save combined spectrum
    np.savetxt(out_file, np.column_stack([wl_combined, flux_combined]),
               fmt="%.6f %.6e",
               header="Wavelength(Angstrom)  Flux(erg/cm2/s/Angstrom)")
    print(f"Combined spectrum saved to {out_file}")

    # Plot 1: Combined spectrum (log scale, 6200–8400 Å)
    mask_combined = (wl_combined >= 6200) & (wl_combined <= 8400)
    plt.figure(figsize=(10, 6))
    plt.plot(wl_combined[mask_combined], flux_combined[mask_combined], label="Combined", color="black", linewidth=1.2)
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux (erg/cm²/s/Å)")
    plt.yscale("log")
    plt.title("Combined Spectrum (Log Scale, 6200–8400 Å)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(combined_plot, dpi=300)
    plt.close()
    print(f"Combined spectrum plot saved to {combined_plot}")

    # Plot 2: H-alpha input spectrum (full range, log y)
    plt.figure(figsize=(10, 6))
    plt.plot(wl2, flux2, color="red", label=f"{base2}")
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux (erg/cm²/s/Å)")
    plt.yscale("log")
    plt.title(f"H-alpha Input Spectrum: {base2}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(halpha_plot, dpi=300)
    plt.close()
    print(f"H-alpha input plot saved to {halpha_plot}")

    # Plot 3: Both inputs together (log y, 5000–8000 Å)
    mask_inputs_cifist = (wl1 >= 5000) & (wl1 <= 8000)
    mask_inputs_halpha = (wl2 >= 5000) & (wl2 <= 8000)
    plt.figure(figsize=(10, 6))
    plt.plot(wl1[mask_inputs_cifist], flux1[mask_inputs_cifist], label=f"{base1}", alpha=0.7)
    plt.plot(wl2[mask_inputs_halpha], flux2[mask_inputs_halpha], label=f"{base2}", alpha=0.7)
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux (erg/cm²/s/Å)")
    plt.yscale("log")
    plt.title(f"Input Spectra (Log Scale, 5000–8000 Å): {base1} & {base2}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(inputs_plot, dpi=300)
    plt.close()
    print(f"Input spectra plot saved to {inputs_plot}")

if __name__ == "__main__":
    main()
