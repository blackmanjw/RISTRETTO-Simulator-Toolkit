import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import CubicSpline
import pandas as pd
import pickle
import glob
from scipy import signal
from barycorrpy import get_BC_vel , exposure_meter_BC_vel
from astropy.time import Time
from scipy.interpolate import interp1d
import plotly.express as px
from expecto import get_spectrum
from PyAstronomy import pyasl
import lmfit
from lmfit import Model
import corner
import random
from datetime import datetime, timedelta
import parameters
from scipy.optimize import least_squares
from scipy.interpolate import interp1d

# ============================================================================
# SPECTRAL EXTRACTION AND MERGING PIPELINE FOR ECHELLE SPECTROGRAPH DATA
# ============================================================================
# This script performs optimal extraction of echelle spectra using pre-computed
# order traces and wavelength solutions. It implements flat-field correction,
# error propagation, and order merging with weighted averaging in overlapping regions.
# ============================================================================

# Base directory containing all necessary calibration files
base_path = '/Users/odin/toolkit/'

# ============================================================================
# LOAD CALIBRATION DATA
# ============================================================================

# Load order traces: polynomial fits describing the location of each spectral order
# on the 2D detector for each fiber
orders = pickle.load(open(base_path + 'pickleFilesMaddalena/orders_from_fit.pickle', "rb"))

# Load wavelength solution: polynomial coefficients mapping pixel position to wavelength
# for each order and fiber
waveSol = pickle.load(open(base_path + 'pickleFilesMaddalena/waveSol_from_HDF.pickle', "rb"))

# Load master flat field: normalized response correction and spatial profiles
# derived from flat-field lamp exposures
flat_pickle = pickle.load(open(base_path + 'masterFlat/master.pickle', "rb"))

# ============================================================================
# EXTRACTION PARAMETERS
# ============================================================================

# Width of the extraction aperture in the cross-dispersion direction (pixels)
# This defines how many pixels perpendicular to the dispersion are summed
width = 10

# CCD readout noise in electrons (used for variance calculation)
readout_noise = 3

# ============================================================================
# OPTIMAL EXTRACTION FUNCTIONS
# ============================================================================

def get_spec_optimal_all(frame, the_fiber, the_order):
    """
    Perform optimal extraction of a spectrum for a single fiber and order.
    
    This implements the Horne (1986) optimal extraction algorithm, which weights
    each pixel by the spatial profile to maximize S/N ratio. The extraction includes:
    - Flat-field correction (deblaze)
    - Proper variance weighting using the spatial profile
    - Error propagation through all steps
    
    Parameters:
    -----------
    frame : 2D array
        The science image (2D CCD frame) to extract spectrum from
    the_fiber : int
        Fiber ID to extract (identifies which fiber on the spectrograph)
    the_order : int
        Spectral order ID to extract (echelle orders numbered by wavelength)
    
    Returns:
    --------
    my_spec_new : 1D array
        Optimally extracted and flat-fielded spectrum (flux in photons)
    my_spec_error : 1D array
        Propagated uncertainties for each spectral pixel
    """
    
    # Create pixel coordinate array spanning the detector width
    x = np.arange(len(frame))
    
    # Get the polynomial trace describing the central position of this order
    # y(x) gives the row position of the order center at each column x
    y = orders[(orders['fiber_ID']==the_fiber) & (orders['order_ID']==the_order)]['poly'].values[0](x)
    
    # Extract a rectangular aperture around the trace for each column
    # spec[i] contains the cross-dispersion slice at column i
    spec = [frame[int(np.round(y[i] - width/2)) : int(np.round(y[i] + width/2)), x[i]] for i in range(len(x))]
    
    # Compute simple sum of flux in each column (for initial variance estimate)
    sum_spec = np.array([np.sum(spec[x[i]]) for i in range(len(x))])
    
    # Load the normalized spatial profile from the flat field
    # This describes how flux is distributed across the aperture
    profile = flat_pickle[(flat_pickle['fiber_ID']==the_fiber) & (flat_pickle['order_ID']==the_order)]['flat_profile'].values[0]
    
    # Calculate variance for each pixel in the aperture
    # Variance = readout_noise^2 + photon_noise (Poisson statistics)
    # Using maximum(0, ...) prevents negative variance from cosmic rays or artifacts
    variance = [readout_noise**2 + np.maximum(0, sum_spec[x[i]]*profile[x[i]]) for i in range(len(x))]
    
    # ---- FIRST ITERATION OF OPTIMAL EXTRACTION ----
    # Weighted sum using profile/variance weights (numerator of Horne equation)
    nominator_s = np.array([np.sum(spec[x[i]] * profile[x[i]] / variance[x[i]]) for i in range(len(x))])
    
    # Sum of squared profile weights (denominator of Horne equation)
    denominator = np.array([np.sum((profile[x[i]])**2 / variance[x[i]]) for i in range(len(x))])
    
    # Optimally extracted spectrum (first pass)
    my_spec = nominator_s / denominator
    
    # ---- SECOND ITERATION WITH IMPROVED VARIANCE ----
    # Recalculate variance using the extracted spectrum for better photon noise estimate
    variance = [readout_noise**2 + np.maximum(0, my_spec[x[i]]*profile[x[i]]) for i in range(len(x))]
    
    # Re-extract with improved variance weights
    nominator_s = np.array([np.sum(spec[x[i]] * profile[x[i]] / variance[x[i]]) for i in range(len(x))])
    nominator_e = 1.0  # Numerator for error (unity for simple variance propagation)
    denominator = np.array([np.sum((profile[x[i]])**2 / variance[x[i]]) for i in range(len(x))])
    
    # Final optimally extracted spectrum
    my_spec = nominator_s / denominator
    
    # Propagated error from the weighted extraction
    my_spec_error = np.sqrt(nominator_e / denominator)
    
    # ---- FLAT-FIELD CORRECTION (DEBLAZE) ----
    # Load the master flat spectrum (blaze function) for this order
    my_flat = flat_pickle[(flat_pickle['fiber_ID']==the_fiber) & (flat_pickle['order_ID']==the_order)]['flat'].values[0]
    
    # Divide by flat to remove instrumental response and blaze function
    my_spec_new = my_spec / my_flat
    
    # ---- ERROR PROPAGATION THROUGH FLAT-FIELDING ----
    # Load uncertainty in the flat field measurement
    my_flat_error = flat_pickle[(flat_pickle['fiber_ID']==the_fiber) & (flat_pickle['order_ID']==the_order)]['flat_error'].values[0]
    
    # Propagate errors: σ(f/g) = (f/g) * sqrt((σf/f)^2 + (σg/g)^2)
    my_spec_error = (my_spec/my_flat) * np.sqrt((my_spec_error/my_spec)**2 + (my_flat_error/my_flat)**2)
    
    return my_spec_new, my_spec_error


def get_spec_optimal_all_NOBLAZE(frame, the_fiber, the_order):
    """
    Optimal extraction WITHOUT flat-field correction (no deblaze).
    
    This function is identical to get_spec_optimal_all() but returns the
    spectrum before flat-field division. Useful for examining raw blaze
    function or for cases where flat-fielding should be done separately.
    
    Parameters and algorithm are identical to get_spec_optimal_all().
    
    Returns:
    --------
    my_spec : 1D array
        Optimally extracted spectrum WITHOUT flat-field correction
    my_spec_error : 1D array
        Propagated uncertainties (before flat-fielding)
    """
    
    x = np.arange(len(frame))
    y = orders[(orders['fiber_ID']==the_fiber) & (orders['order_ID']==the_order)]['poly'].values[0](x)
    spec = [frame[int(np.round(y[i] - width/2)) : int(np.round(y[i] + width/2)), x[i]] for i in range(len(x))]
    sum_spec = np.array([np.sum(spec[x[i]]) for i in range(len(x))])
    profile = flat_pickle[(flat_pickle['fiber_ID']==the_fiber) & (flat_pickle['order_ID']==the_order)]['flat_profile'].values[0]
    variance = [readout_noise**2 + np.maximum(0, sum_spec[x[i]]*profile[x[i]]) for i in range(len(x))]
    
    nominator_s = np.array([np.sum(spec[x[i]] * profile[x[i]] / variance[x[i]]) for i in range(len(x))])
    denominator = np.array([np.sum((profile[x[i]])**2 / variance[x[i]]) for i in range(len(x))])
    my_spec = nominator_s / denominator
    
    # Second iteration
    variance = [readout_noise**2 + np.maximum(0, my_spec[x[i]]*profile[x[i]]) for i in range(len(x))]
    nominator_s = np.array([np.sum(spec[x[i]] * profile[x[i]] / variance[x[i]]) for i in range(len(x))])
    nominator_e = 1.0
    denominator = np.array([np.sum((profile[x[i]])**2 / variance[x[i]]) for i in range(len(x))])
    my_spec = nominator_s / denominator
    my_spec_error = np.sqrt(nominator_e / denominator)
    
    return my_spec, my_spec_error


def get_wave(waveSol, the_fiber, the_order):
    """
    Retrieve the wavelength solution for a given fiber and order.
    
    The wavelength solution maps pixel position to wavelength using a
    polynomial fit (typically from ThAr lamp calibration).
    
    Parameters:
    -----------
    waveSol : DataFrame
        Wavelength solution lookup table with polynomials
    the_fiber : int
        Fiber ID
    the_order : int
        Spectral order ID
    
    Returns:
    --------
    wavelength : 1D array
        Wavelength in Angstroms for each pixel (0-4095)
    """
    
    # Evaluate the wavelength polynomial at each pixel position
    wavelength = waveSol[(waveSol['fiber_ID']==the_fiber) & (waveSol['order_ID']==the_order)]['poly'].values[0](np.arange(4096))
    
    # Convert to Angstroms (input is in microns or nm)
    return wavelength * 1e4


def blaze_func(flux1, flux_flat, error1, error_flat):
    """
    Apply blaze function correction (flat-fielding) with error propagation.
    
    Note: This function is defined but not used in the main pipeline
    (flat-fielding is done inside get_spec_optimal_all instead).
    
    Parameters:
    -----------
    flux1 : array
        Raw spectrum flux
    flux_flat : array
        Flat field (blaze function)
    error1 : array
        Uncertainties in raw spectrum
    error_flat : array
        Uncertainties in flat field
    
    Returns:
    --------
    flux : array
        Flat-fielded spectrum
    error : array
        Propagated errors
    """
    
    flux = flux1 / flux_flat
    error = np.abs(flux) * np.sqrt((error1/flux1)**2 + (error_flat/flux_flat)**2)
    return flux, error


# ============================================================================
# UTILITY FUNCTIONS FOR ORDER MERGING
# ============================================================================

def elinan(sp1, wavelengths1):
    """
    Remove NaN values from spectrum and wavelength arrays.
    
    NaN values can arise from bad pixels, cosmic rays, or edge effects.
    This function identifies and removes all indices where NaN appears in
    either the flux, error, or wavelength arrays.
    
    Parameters:
    -----------
    sp1 : tuple of 2 arrays
        (flux_array, error_array)
    wavelengths1 : array
        Wavelength array
    
    Returns:
    --------
    sp1 : tuple of 2 arrays
        Cleaned (flux, error) with NaNs removed
    wavelengths1 : array
        Cleaned wavelength array
    """
    
    sp1 = list(sp1)
    
    # Find indices of NaN in flux array
    a = np.where(np.isnan(sp1[0]))[0]
    
    # Find indices of NaN in error array
    b = np.where(np.isnan(sp1[1]))[0]
    
    # Find indices of NaN in wavelength array
    c = np.where(np.isnan(wavelengths1))[0]
    
    # Combine all NaN indices (unique values only)
    d = np.unique(np.concatenate([a, b, c]))
    
    # Remove NaN indices from all arrays
    sp1[0] = np.delete(sp1[0], d)
    sp1[1] = np.delete(sp1[1], d)
    wavelengths1 = np.delete(wavelengths1, d)
    
    return sp1, wavelengths1


def weighted_avg_and_std(fluxes, weights):
    """
    Calculate weighted average and uncertainty for combining overlapping orders.
    
    When two spectral orders overlap, we combine them using inverse-variance
    weighting to maximize S/N. This function computes the weighted mean and
    its uncertainty.
    
    Parameters:
    -----------
    fluxes : array
        Array of flux values to combine (typically 2 values from 2 orders)
    weights : array
        Weight for each flux value (typically 1/sigma^2)
    
    Returns:
    --------
    weighted_avg : float
        Weighted mean flux
    weighted_error : float
        Uncertainty in weighted mean (sqrt of weighted variance)
    """
    
    # Weighted average: sum(w_i * f_i) / sum(w_i)
    weighted_avg = np.sum(weights * fluxes) / np.sum(weights)
    
    # Variance of weighted average: 1 / sum(w_i)
    weighted_variance = 1 / np.sum(weights)
    
    return weighted_avg, np.sqrt(weighted_variance)


# ============================================================================
# ORDER MERGING FUNCTION
# ============================================================================

def merge_all_orders(frame, the_fiber):
    """
    Merge all spectral orders into a single continuous spectrum.
    
    Echelle spectrographs produce multiple overlapping spectral orders.
    This function:
    1. Extracts each order individually
    2. Interpolates all orders onto a common wavelength grid (log-spaced)
    3. Combines overlapping regions using inverse-variance weighted averaging
    4. Stitches orders together to create a continuous spectrum
    
    The wavelength grid is logarithmically spaced to provide constant velocity
    resolution, which is optimal for radial velocity measurements.
    
    Parameters:
    -----------
    frame : 2D array
        Science image to extract and merge
    the_fiber : int
        Fiber ID to extract
    
    Returns:
    --------
    wave_merged : 1D array
        Wavelength array for merged spectrum (Angstroms)
    spectrum_merged : tuple of 2 arrays
        (merged_flux, merged_error) for the continuous spectrum
    """
    
    # Extraction aperture width (must match global parameter)
    width = 10
    
    # Initialize merged spectrum placeholders
    wave_merged = None
    spectrum_merged = None
    
    # ---- DEFINE LOGARITHMIC WAVELENGTH GRID ----
    # Sampling in velocity space (km/s per pixel)
    sampling = 0.86
    
    # Create log-spaced wavelength grid from 610 nm to 845 nm
    # Spacing: Δλ/λ = v/c where v is the velocity sampling
    logwave = np.arange(np.log10(610), np.log10(845), np.log10(1. + sampling/2.99792458e5))
    
    # Convert back to linear wavelength scale
    wavelength_grid_angstroms = 10**logwave
    
    # Trim edges to avoid interpolation artifacts (in Angstroms)
    wave = wavelength_grid_angstroms[100:-100] * 10
    
    # ---- ITERATE THROUGH ORDERS FROM BLUE TO RED ----
    # Orders are numbered 123 (bluest) to 90 (reddest)
    start_order = 123
    
    for i in range(start_order, 90, -1):
        
        # Extract spectrum for current order
        spec2 = get_spec_optimal_all(frame, the_fiber, i)
        wave2 = get_wave(waveSol, the_fiber, i)
        
        # ---- INTERPOLATE TO COMMON WAVELENGTH GRID ----
        # Find wavelength range covered by this order in the common grid
        overlapped = np.where((wave >= wave2[0]) & (wave <= wave2[-1]))[0]
        
        # Cubic spline interpolation of flux and error onto common grid
        spec2 = (CubicSpline(wave2, spec2[0])(wave[overlapped]), 
                 CubicSpline(wave2, spec2[1])(wave[overlapped]))
        wave2 = wave[overlapped]
        
        # Extract flux and error for current order
        flux2, error2 = spec2[0], spec2[1]
        
        # ---- HANDLE FIRST ITERATION (BLUEST ORDER) ----
        if i == start_order:
            # For the first order, also extract the next bluer order (124)
            spec1 = get_spec_optimal_all(frame, the_fiber, start_order + 1)
            wave1 = get_wave(waveSol, the_fiber, start_order + 1)
            
            # Interpolate to common grid
            overlapped = np.where((wave >= wave1[0]) & (wave <= wave1[-1]))[0]
            spec1 = (CubicSpline(wave1, spec1[0])(wave[overlapped]),
                     CubicSpline(wave1, spec1[1])(wave[overlapped]))
            wave1 = wave[overlapped]
        else:
            # For subsequent iterations, use the previously merged spectrum
            spec1 = spectrum_merged
            wave1 = wave_merged
        
        # Extract flux and error for previous/merged spectrum
        flux1, error1 = spec1[0], spec1[1]
        
        # ---- IDENTIFY OVERLAPPING REGION ----
        # Find wavelength range where both orders have data
        overlap_start = np.max([wave1[0], wave2[0]])
        overlap_end = np.min([wave1[-1], wave2[-1]])
        
        # Extract common wavelength points in overlap region
        common_wavelength = wave[np.where((wave >= overlap_start) & (wave <= overlap_end))[0]]
        
        # ---- INTERPOLATE BOTH SPECTRA IN OVERLAP REGION ----
        # Create splines for flux and error in overlapping region
        over1 = CubicSpline(wave1, flux1)
        over_err1 = CubicSpline(wave1, error1)
        over2 = CubicSpline(wave2, flux2)
        over_err2 = CubicSpline(wave2, error2)
        
        # ---- CALCULATE INVERSE-VARIANCE WEIGHTS ----
        # Weight each point by 1/σ^2 to minimize combined uncertainty
        weights1 = 1.0 / (over_err1(common_wavelength)**2)
        weights2 = 1.0 / (over_err2(common_wavelength)**2)
        
        # ---- MERGE OVERLAPPING REGION ----
        # For each wavelength point, compute weighted average of both orders
        merged_flux = np.array([
            weighted_avg_and_std(
                np.array([over1(common_wavelength[k]), over2(common_wavelength[k])]),
                np.array([weights1[k], weights2[k]])
            )[0] for k in range(len(weights1))
        ])
        
        # Propagate errors through weighted averaging
        merged_error = np.array([
            weighted_avg_and_std(
                np.array([over1(common_wavelength[k]), over2(common_wavelength[k])]),
                np.array([weights1[k], weights2[k]])
            )[1] for k in range(len(weights1))
        ])
        
        # ---- CONCATENATE MERGED SPECTRUM ----
        # Combine: [order1 non-overlap] + [merged overlap] + [order2 non-overlap]
        wave_merged = np.concatenate((
            wave1[wave1 < overlap_start],
            common_wavelength,
            wave2[wave2 > overlap_end]
        ))
        
        spectrum_merged = (
            np.concatenate((
                flux1[wave1 < overlap_start],
                merged_flux,
                flux2[wave2 > overlap_end]
            )),
            np.concatenate((
                error1[wave1 < overlap_start],
                merged_error,
                error2[wave2 > overlap_end]
            ))
        )
    
    return wave_merged, spectrum_merged


# ============================================================================
# MAIN EXECUTION: EXTRACT AND PLOT SPECTRUM
# ============================================================================

def main():
    """
    Main function to run spectral extraction from command line.
    
    Usage:
        python analysis.py <filename> [fiber_number]
    
    Arguments:
        filename: Name of FITS file (without .fits extension) in pychelle_output/
        fiber_number: Optional fiber ID to extract (default: 1)
    
    Example:
        python analysis.py star_3600_s 1
        python analysis.py star_3600_s     # defaults to fiber 1
    """
    import sys
    import os
    
    # ---- PARSE COMMAND LINE ARGUMENTS ----
    if len(sys.argv) < 2:
        print("Usage: python analysis.py <filename> [fiber_number]")
        print("Example: python analysis.py star_3600_s 1")
        sys.exit(1)
    
    # Get filename from first argument (without .fits extension)
    filename = sys.argv[1]
    
    # Get fiber number from second argument, default to 1 if not provided
    the_fiber = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    
    # ---- CREATE ANALYSIS DIRECTORY IF IT DOESN'T EXIST ----
    analysis_dir = base_path + 'analysis'
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
        print(f"Created directory: {analysis_dir}")
    
    # Construct full path to FITS file
    file_fits = base_path + f'pychelle_output/{filename}.fits'
    
    print(f"Processing file: {file_fits}")
    print(f"Extracting fiber: {the_fiber}")
    
    # Read 2D image and subtract bias level (250 ADU)
    try:
        spectrum_2D = (fits.getdata(file_fits).astype('float') - 250)
    except FileNotFoundError:
        print(f"Error: File not found: {file_fits}")
        sys.exit(1)
    
    # Extract and merge all orders for specified fiber
    print("Merging spectral orders...")
    wave_merged, merged_spec = merge_all_orders(spectrum_2D, the_fiber)
    
    # ---- PLOT MERGED SPECTRUM ----
    plt.figure(figsize=(10, 6))
    plt.title(f'Merged Spectrum - {filename} (Fiber {the_fiber})')
    plt.plot(wave_merged, merged_spec[0])
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Flux (photons)')
    plt.xlim(6200, 6500)  # Show small wavelength region
    plt.grid(alpha=0.3)
    plt.tight_layout()
    output_file = os.path.join(analysis_dir, f'{filename}_fiber{the_fiber}_merged.png')
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.show()
    
    # ---- EXTRACT AND PLOT SINGLE ORDER ----
    the_order = 117  # Select spectral order as example
    
    print(f"Extracting single order {the_order}...")
    # Extract single order without merging
    spectrum1D_order = get_spec_optimal_all(spectrum_2D, the_fiber, the_order)
    wavelengths = get_wave(waveSol, the_fiber, the_order)
    
    # Plot single order spectrum
    plt.figure(figsize=(10, 6))
    plt.title(f'Spectrum for Fiber {the_fiber}, Order {the_order} - {filename}')
    plt.plot(wavelengths, spectrum1D_order[0])
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Flux (photons)')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    output_file = os.path.join(analysis_dir, f'{filename}_fiber{the_fiber}_order{the_order}.png')
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=150)
    print(f"Saved: {output_file}")
    plt.show()
    
    print("Processing complete!")


# Run main function when script is executed
if __name__ == "__main__":
    main()
