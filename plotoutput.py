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

    plt.figure(figsize=(12, 6))

    for f in files:
        try:
            # Skip header row, assume two numeric columns
            wl, flux = np.loadtxt(f, delimiter=",", unpack=True, skiprows=1)
            label = os.path.basename(f).replace("_fiber1.csv", "")
            plt.plot(wl, flux, label=label, alpha=0.8)
        except Exception as e:
            print(f"Skipping {f}: {e}")

    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Flux")
    plt.title("All fiber1 spectra")

    # Place legend outside plot on the right
    plt.legend(
        fontsize="small",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0
    )

    plt.yscale("log")   # remove this if you want linear flux
    plt.tight_layout(rect=[0, 0, 0.8, 1])  # leave space on right for legend

    # Save plot as PNG inside output folder
    save_path = os.path.join(output_dir, "plotoutput.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved to {save_path}")

    plt.show()

if __name__ == "__main__":
    main()
