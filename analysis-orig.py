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

# Extraction spectrum Pyechelle

base_path = '/Users/odin/toolkit'

#################### Optimal extraction ####################

orders=pickle.load(open(base_path + 'pickleFilesMaddalena/orders_from_fit.pickle', "rb"))
waveSol= pickle.load(open(base_path + 'pickleFilesMaddalena/waveSol_from_HDF.pickle', "rb"))
flat_pickle=pickle.load(open(base_path + 'masterFlat/master.pickle', "rb"))

################################# Create profile from master flat field ######################################

width=10   # width of the spectrum in the cross-dispersion direction in pixels
readout_noise = 3

#################### Extract the spectrum of a single fiber and a single order ####################

def get_spec_optimal_all(frame, the_fiber, the_order):
    x                = np.arange(len(frame))
    y                = orders[(orders['fiber_ID']==the_fiber) & (orders['order_ID']==the_order)]['poly'].values[0](x)                          
    spec             = [frame[int(np.round(y[i] - width/2)) : int(np.round(y[i] + width/2)), x[i]] for i in range(len(x))]
    sum_spec         =  np.array([np.sum(spec[x[i]]) for i in range(len(x))])
    profile          = flat_pickle[ (flat_pickle['fiber_ID']==the_fiber) & (flat_pickle['order_ID']==the_order) ]['flat_profile'].values[0]   
    variance         = [readout_noise**2 + np.maximum(0,sum_spec[x[i]]*profile[x[i]]) for i in range(len(x))]  #[readout_noise**2 + sum_spec[x[i]]*profile[x[i]] for i in range(len(x))] #

    #variance         = [readout_noise**2 + np.maximum(0,frame[int(np.round(y[i] - width/2)) : int(np.round(y[i] + width/2)), x[i]])*np.array(profile)[x[i]] for i in range(len(x))]

    nominator_s      = np.array([ np.sum( spec[x[i]] * profile[x[i]] / variance[x[i]] ) for i in range(len(x))])
    denominator      = np.array([ np.sum( (profile[x[i]])**2 / variance[x[i]] ) for i in range(len(x))])
    my_spec          = nominator_s / denominator

    variance         = [readout_noise**2 + np.maximum(0,my_spec[x[i]]*profile[x[i]]) for i in range(len(x))]  
    nominator_s      = np.array([ np.sum( spec[x[i]] * profile[x[i]] / variance[x[i]] ) for i in range(len(x))])
    nominator_e      = 1. 
    denominator      = np.array([ np.sum( (profile[x[i]])**2 / variance[x[i]] ) for i in range(len(x))])
    my_spec          = nominator_s / denominator
    my_spec_error    = np.sqrt(nominator_e / denominator)
    
    my_flat=flat_pickle[(flat_pickle['fiber_ID']==the_fiber) & (flat_pickle['order_ID']==the_order)]['flat'].values[0]
    my_spec_new=my_spec/my_flat
    ################################################ Error propagation of the flat field ##########################################
    my_flat_error=flat_pickle[(flat_pickle['fiber_ID']==the_fiber) & (flat_pickle['order_ID']==the_order)]['flat_error'].values[0]
    my_spec_error=(my_spec/my_flat)*np.sqrt((my_spec_error/my_spec)**2+(my_flat_error/my_flat)**2)
    ###############################################################################################################################
    return my_spec_new, my_spec_error#/my_flat

################### Not deblazed spectrum ####################

def get_spec_optimal_all_NOBLAZE(frame, the_fiber, the_order):
    x                = np.arange(len(frame))
    y                = orders[(orders['fiber_ID']==the_fiber) & (orders['order_ID']==the_order)]['poly'].values[0](x)                          
    spec             = [frame[int(np.round(y[i] - width/2)) : int(np.round(y[i] + width/2)), x[i]] for i in range(len(x))]
    sum_spec         =  np.array([np.sum(spec[x[i]]) for i in range(len(x))])
    profile          = flat_pickle[ (flat_pickle['fiber_ID']==the_fiber) & (flat_pickle['order_ID']==the_order) ]['flat_profile'].values[0]   
    variance         = [readout_noise**2 + np.maximum(0,sum_spec[x[i]]*profile[x[i]]) for i in range(len(x))]  #[readout_noise**2 + sum_spec[x[i]]*profile[x[i]] for i in range(len(x))] #

    #variance         = [readout_noise**2 + np.maximum(0,frame[int(np.round(y[i] - width/2)) : int(np.round(y[i] + width/2)), x[i]])*np.array(profile)[x[i]] for i in range(len(x))]

    nominator_s      = np.array([ np.sum( spec[x[i]] * profile[x[i]] / variance[x[i]] ) for i in range(len(x))])
    denominator      = np.array([ np.sum( (profile[x[i]])**2 / variance[x[i]] ) for i in range(len(x))])
    my_spec          = nominator_s / denominator

    variance         = [readout_noise**2 + np.maximum(0,my_spec[x[i]]*profile[x[i]]) for i in range(len(x))]  
    nominator_s      = np.array([ np.sum( spec[x[i]] * profile[x[i]] / variance[x[i]] ) for i in range(len(x))])
    nominator_e      = 1. 
    denominator      = np.array([ np.sum( (profile[x[i]])**2 / variance[x[i]] ) for i in range(len(x))])
    my_spec          = nominator_s / denominator
    my_spec_error    = np.sqrt(nominator_e / denominator)

    return my_spec, my_spec_error

#################### Extract the wavelength solution of a given order and fiber ####################

def get_wave(waveSol, the_fiber, the_order):
    wavelength       = waveSol[(waveSol['fiber_ID']==the_fiber) & (waveSol['order_ID']==the_order)]['poly'].values[0](np.arange(4096))
    return wavelength*1e4 # in Angstroms

def blaze_func(flux1,flux_flat,error1,error_flat):
    flux=flux1/flux_flat
    error=np.abs(flux)*np.sqrt((error1/flux1)**2+(error_flat/flux_flat)**2)
    return flux,error

#################### Useful functions for order merging ####################

# Eliminate NaN values from the spectrum and the wavelength
def elinan(sp1,wavelengths1):
    sp1=list(sp1)
    a=np.where(np.isnan(sp1[0]))[0]
    b=np.where(np.isnan(sp1[1]))[0]
    c=np.where(np.isnan(wavelengths1))[0]
    d=np.unique(np.concatenate([a,b,c]))
    sp1[0]=np.delete(sp1[0],d)
    sp1[1]=np.delete(sp1[1],d)
    wavelengths1=np.delete(wavelengths1,d)
    return sp1,wavelengths1

def weighted_avg_and_std(fluxes, weights):
    weighted_avg = np.sum(weights * fluxes) / np.sum(weights)
    weighted_variance=1/np.sum(weights)
    return weighted_avg, np.sqrt(weighted_variance)

#################### Useful functions for order merging ####################

# Eliminate NaN values from the spectrum and the wavelength
def elinan(sp1,wavelengths1):
    sp1=list(sp1)
    a=np.where(np.isnan(sp1[0]))[0]
    b=np.where(np.isnan(sp1[1]))[0]
    c=np.where(np.isnan(wavelengths1))[0]
    d=np.unique(np.concatenate([a,b,c]))
    sp1[0]=np.delete(sp1[0],d)
    sp1[1]=np.delete(sp1[1],d)
    wavelengths1=np.delete(wavelengths1,d)
    return sp1,wavelengths1

def weighted_avg_and_std(fluxes, weights):
    weighted_avg = np.sum(weights * fluxes) / np.sum(weights)
    weighted_variance=1/np.sum(weights)
    return weighted_avg, np.sqrt(weighted_variance)

#################### Function to merge all orders ####################

def merged_spectrum(frame, the_fiber):
    width=10
    # Initialize variables for the merged spectrum
    wave_merged = None
    spectrum_merged = None

    sampling=0.86		# km/s/pixel
    logwave=np.arange(np.log10(610),np.log10(845),np.log10(1.+sampling/2.99792458e5))
    wavelength_grid_angstroms=10**logwave
    wave=wavelength_grid_angstroms[100:-100]*10 # in Angstroms

    # Iterate over wavelength orders from 123 to 90
    start_order=123
    for i in range(start_order, 90, -1): #90
        # Get spectrum and wavelength for the current order
        spec2 = get_spec_optimal_all(frame, the_fiber,i)
        wave2 = get_wave(waveSol, the_fiber, i)
        # Cubicspline interpolation with a wavelength grid constant in radial velocities (interpolation of both the flux and the error)
        overlapped = np.where((wave >= wave2[0]) & (wave <= wave2[-1]))[0]
        spec2=CubicSpline(wave2,spec2[0])(wave[overlapped]), CubicSpline(wave2,spec2[1])(wave[overlapped]) 
        wave2=wave[overlapped]

        # Extract flux and error from the current spectrum
        flux2, error2 = spec2[0], spec2[1]

        # For the first iteration (i=123), obtain the spectrum and wavelength for the next order (124)
        if i == start_order:
            spec1 = get_spec_optimal_all(frame, the_fiber, start_order+1)
            wave1 = get_wave(waveSol, the_fiber, start_order+1)
            # Cubicspline interpolation with a wavelength grid constant in radial velocities (interpolation of both the flux and the error)
            overlapped = np.where((wave >= wave1[0]) & (wave <= wave1[-1]))[0]
            spec1=CubicSpline(wave1,spec1[0])(wave[overlapped]), CubicSpline(wave1,spec1[1])(wave[overlapped])
            wave1=wave[overlapped]
        else:
            # For subsequent iterations, use the previously merged spectrum and wavelength
            spec1 = spectrum_merged
            wave1 = wave_merged

        # Extract flux and error from the previous spectrum
        flux1, error1 = spec1[0], spec1[1]

        # Identify the overlapping region between the two wavelength orders
        overlap_start = np.max([wave1[0], wave2[0]])
        overlap_end = np.min([wave1[-1], wave2[-1]])

        common_wavelength=wave[np.where((wave>=overlap_start) & (wave<=overlap_end))[0]]

        # Create cubic splines for the overlapping regions of both spectra
        over1 = CubicSpline(wave1, flux1)
        over_err1 = CubicSpline(wave1, error1)
        over2 = CubicSpline(wave2, flux2)
        over_err2 = CubicSpline(wave2, error2)

        # Calculate weights based on the inverse of errors for both spectra
        weights1 = 1.0 / (over_err1(common_wavelength)**2)
        weights2 = 1.0 / (over_err2(common_wavelength)**2)

        # Calculate the merged flux and error using weighted averages and standard deviations
        merged_flux = np.array([weighted_avg_and_std(np.array([over1(common_wavelength[k]), over2(common_wavelength[k])]), np.array([weights1[k], weights2[k]]))[0] for k in range(len(weights1))])
        merged_error = np.array([weighted_avg_and_std(np.array([over1(common_wavelength[k]), over2(common_wavelength[k])]), np.array([weights1[k], weights2[k]]))[1] for k in range(len(weights1))])

        # Concatenate the wavelength, flux, and error arrays to form the merged spectrum
        wave_merged = np.concatenate((wave1[wave1 < overlap_start], common_wavelength, wave2[wave2 > overlap_end]))

        spectrum_merged = (
            np.concatenate((flux1[wave1 < overlap_start], merged_flux, flux2[wave2 > overlap_end])),
            np.concatenate((error1[wave1 < overlap_start], merged_error, error2[wave2 > overlap_end]))
        )
        
    return wave_merged,spectrum_merged

file_fits = base_path + 'outPyechelle/star_3600_s.fits'
# Open merged file
spectrum_2D=(fits.getdata(file_fits).astype('float')-250)
wave_merged, merged_spectrum = merged_spectrum(spectrum_2D, 1)

plt.figure(figsize=(10, 6))
plt.title('Merged Spectrum')
plt.plot(wave_merged, merged_spectrum[0])
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (photons)')
plt.xlim(6200,6500)
plt.legend()
plt.show()

# Extract 1 order
the_order=117
the_fiber=1
spectrum1D_order=get_spec_optimal_all(spectrum_2D,the_fiber, the_order)
wavelengths=get_wave(waveSol, the_fiber, the_order)
# plot
plt.figure(figsize=(10, 6))
plt.title(f'Spectrum for Fiber {the_fiber}, Order {the_order}')
plt.plot(wavelengths, spectrum1D_order[0])
plt.xlabel('Wavelength (Angstroms)')
plt.ylabel('Flux (photons)')
plt.legend()
plt.show()
