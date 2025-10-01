#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

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
    plt.legend(fontsize="small", ncol=2)
    plt.yscale("log")   # remove this if you want linear flux
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
