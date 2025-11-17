import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from matplotlib.ticker import ScalarFormatter
import os
import argparse
import smplotlib

def main():
    parser = argparse.ArgumentParser(description="Resample spectrum to a target resolution.")
    parser.add_argument("filename", help="Name of the spectrum file in the input folder")
    parser.add_argument("--R", type=int, default=140000, help="Target resolution (default: 140000)")
    parser.add_argument("--no-show", action="store_true", help="Do not display the plot interactively")
    args = parser.parse_args()

    input_folder = "input"
    output_folder = "input"

    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    path = os.path.join(input_folder, args.filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"File '{path}' not found.")

    # Load spectrum
    wavelength, flux = np.loadtxt(path, unpack=True)

    # Compute sigma for Gaussian smoothing
    delta_lambda = wavelength / args.R
    sigma_pixels = delta_lambda / (wavelength[1] - wavelength[0]) / 2.355
    sigma_median = np.median(sigma_pixels)

    # Smooth flux
    flux_smoothed = gaussian_filter1d(flux, sigma=sigma_median)

    # Plot result
    plt.figure(figsize=(10, 5))
    plt.plot(wavelength, flux, label="Original")
    plt.plot(wavelength, flux_smoothed, label=f"Resampled (R={args.R})")
    plt.legend()
    plt.xlim(6200, 8400)

    # Disable scientific notation on axes
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=False))
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.xaxis.get_major_formatter().set_useOffset(False)

    # Conditional y-axis limits
    if "star" in args.filename.lower():
        plt.ylim(bottom=100000)
    else:
        plt.ylim(0.00001, 3e3)

    plt.yscale('log')
    plt.xlabel(r"Wavelength ($\mathrm{\AA}$)")
    plt.ylabel(r"Flux (erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}^{-1}$)")

    # Determine if _orig should be included
    special_objects = ["2MJ1612b", "PDS70b", "PDS70c", "WISPIT2b"]
    use_orig = not any(obj in args.filename for obj in special_objects)

    # Save figure
    name, ext = os.path.splitext(args.filename)
    if use_orig:
        plot_filename = os.path.join(output_folder, f"{name}_{args.R}_orig.png")
    else:
        plot_filename = os.path.join(output_folder, f"{name}_{args.R}.png")
    plt.savefig(plot_filename, bbox_inches='tight', dpi=300)

    # Show plot only if not suppressed
    if not args.no_show:
        plt.show()

    # Save resampled spectrum
    if use_orig:
        output_filename = os.path.join(output_folder, f"{name}_{args.R}_orig{ext}")
    else:
        output_filename = os.path.join(output_folder, f"{name}_{args.R}{ext}")
    np.savetxt(
        output_filename,
        np.column_stack([wavelength, flux_smoothed]),
        fmt=["%.3f", "%.6e"],  # wavelength with 3 decimals, flux in scientific notation
        header=f"Resampled spectrum at R={args.R}"
    )
    print(f"Saved resampled spectrum to {output_filename}")
    print(f"Saved plot to {plot_filename}")

if __name__ == "__main__":
    main()
