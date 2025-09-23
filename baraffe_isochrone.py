import numpy as np
from scipy.interpolate import RegularGridInterpolator

# Conversion constants
MSUN_TO_MJUP = 1047.3486


def parse_track_file(filename):
    """
    Parse stellar evolution track file into structured arrays.
    Returns:
        ages (array, Myr), masses (array, Mjup), data_dict (dict of 2D arrays)
    """
    ages = []
    mass_grid = []
    data_blocks = {}

    current_age = None
    columns = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # Detect new age block
            if line.startswith("t (Gyr)"):
                current_age = float(line.split('=')[1]) * 1000.0  # convert Gyr -> Myr
                ages.append(current_age)
                mass_grid.append([])
                continue

            # Skip header lines
            if line.startswith("M/Ms"):
                columns = line.split()
                # Initialize storage for each column
                for col in columns:
                    if col not in data_blocks:
                        data_blocks[col] = []
                continue

            # Parse numerical row
            parts = line.split()
            if len(parts) == len(columns):
                row = {col: float(val) for col, val in zip(columns, parts)}
                mass_val = row["M/Ms"] * MSUN_TO_MJUP  # convert Msun -> Mjup
                mass_grid[-1].append(mass_val)
                for col, val in row.items():
                    data_blocks[col].append(val)

    # Convert to arrays
    ages = np.array(ages)
    unique_masses = np.unique(np.concatenate(mass_grid))

    # Build 2D grids for each property
    n_age = len(ages)
    n_mass = len(unique_masses)
    data_dict = {}

    for col in columns:
        arr = np.full((n_age, n_mass), np.nan)
        idx = 0
        for i, age in enumerate(ages):
            for j, mass in enumerate(mass_grid[i]):
                arr[i, j] = data_blocks[col][idx]
                idx += 1
        data_dict[col] = arr

    return ages, unique_masses, data_dict


def build_interpolators(ages, masses, data_dict):
    """
    Build interpolators for all columns.
    Returns dict of interpolators keyed by column name.
    """
    interpolators = {}
    for col, grid in data_dict.items():
        interpolators[col] = RegularGridInterpolator(
            (ages, masses), grid, bounds_error=False, fill_value=None
        )
    return interpolators


def get_params(interpolators, age, mass, columns=("Teff", "g")):
    """
    Query interpolators for given age (Myr) and mass (Mjup).
    Returns dict with requested columns.
    """
    point = np.array([age, mass])
    return {col: float(interpolators[col](point)) for col in columns}


# Example usage:
# ages, masses, data_dict = parse_track_file("track.dat")
# interps = build_interpolators(ages, masses, data_dict)
# result = get_params(interps, 5.1, 5.3)  # one line query
# print(result)  # {"Teff": ..., "g": ...}
