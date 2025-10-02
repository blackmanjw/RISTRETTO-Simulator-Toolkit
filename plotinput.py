#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import smplotlib  # keep this since you added it

def main():
    input_dir = "input"

    # Original pattern
    pattern1 = os.path.join(input_dir, "*140000_Ha_60_12.txt")
    files1 = glob.glob(pattern1)

    # New pattern: *140000.txt containing 'star'
    pattern2 = os.path.join(input_dir, "*140000.txt")
    files2 = [f for f in glob.glob(pattern2) if "star" in os.path.basename(f)]

    # Combine both lists and sort
    files = sorted(files1 + files2)

    if not files:
        print("No matching *.txt files found in 'input/'")
        return

    fig, ax = plt.subplots(figsize=(12, 6))

    for f in files:
        try:
            # Load while ignoring '#' comment lines
            wl, flux = np.loadtxt(f, unpack=True, comments="#")
            # Convert Å → nm
            wl = wl / 10.0

            # Mask: only keep 600–850 nm
            mask = (wl >= 600.0) & (wl <= 850.0)
            wl = wl[mask]
            flux = flux[mask]

            label = os.path.basename(f).replace(".txt", "")
            ax.plot(wl, flux, label=label, alpha=0.8)
        except Exception as e:
            print(f"Skipping {f}: {e}")

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Flux (erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)")
    ax.set_title("All input spectra (600-850 nm)")
    ax.set_yscale("log")

    # Place legend under plot with multiple rows
    ax.legend(
        fontsize="small",
        loc="upper center",
        bbox_to_anchor=(0.5, -0.3),
        ncol=2,
        frameon=False
    )

    # Leave extra space below for legend
    fig.subplots_adjust(bottom=0.35)

    # Save plot as PNG inside input folder
    save_path = os.path.join(input_dir, "plotinputs.png")
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved to {save_path}")

    plt.show()

if __name__ == "__main__":
    main()

