#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import smplotlib  # since you mentioned adding this

def main():
    output_dir = "output"
    pattern = os.path.join(output_dir, "*fiber1.csv")
    files = sorted(glob.glob(pattern))

    if not files:
        print("No *fiber1.csv files found in 'output/'")
        return

    fig, ax = plt.subplots(figsize=(12, 6))

    for f in files:
        try:
            # Skip header row, assume two numeric columns
            wl, flux = np.loadtxt(f, delimiter=",", unpack=True, skiprows=1)
            label = os.path.basename(f).replace("_fiber1.csv", "")
            ax.plot(wl, flux, label=label, alpha=0.8)
        except Exception as e:
            print(f"Skipping {f}: {e}")

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Flux (photons/s)")
    ax.set_title("All spectra for Pyechelle input")
    ax.set_yscale("log")   # remove this if you want linear flux

    # Place legend under plot with more rows (fewer columns)
    legend = ax.legend(
        fontsize="small",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.3),  # push further down
        ncol=2,                      # fewer columns → more rows
        frameon=False
    )

    # Leave space for legend below
    fig.subplots_adjust(bottom=0.35)

    # Save plot as PNG inside output folder
    save_path = os.path.join(output_dir, "plotoutput.png")
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved to {save_path}")

    plt.show()

if __name__ == "__main__":
    main()
