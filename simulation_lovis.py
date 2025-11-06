
import os,sys,string,glob
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import smplotlib
import os

filelist=glob.glob('harps_pds70/*.fits')
filelist.sort()

spe=np.zeros([len(filelist),12000])*0.
exptime=np.zeros(len(filelist))*0.
seeing=np.zeros(len(filelist))*0.
airmass=np.zeros(len(filelist))*0.
SNR=np.zeros(len(filelist))*0.

for i in range(len(filelist)):
	f=fits.open(filelist[i])
	data=f[1].data[0]
	spe[i]=data[1][(data[0]>6500.)&(data[0]<=6620.)]
	exptime[i]=f[0].header['EXPTIME']
	seeing[i]=(f[0].header['HIERARCH ESO TEL AMBI FWHM START']+f[0].header['HIERARCH ESO TEL AMBI FWHM END'])/2.
	airmass[i]=(f[0].header['HIERARCH ESO TEL AIRM START']+f[0].header['HIERARCH ESO TEL AIRM END'])/2.
	SNR[i]=f[0].header['HIERARCH ESO DRS SPE EXT SN67']

wave=data[0][(data[0]>6500.)&(data[0]<=6620.)]
# SNR=spe[:,(wave>6583.)&(wave<6585.)].mean(1)/spe[:,(wave>6583.)&(wave<6585.)].std(1)
rv=(wave-6562.82)/6562.82*3e5

# Flux calibration
# From flux per unit time vs seeing plot:
# Measured mean flux in [6583-6585] A = 0.8 e-/s/bin at an estimated instrument throughput of 7%
tel_surface=88564.
flux_density=0.8*(6.626e-34*2.99792458e8/6584.*1e10)*1e7/0.01/tel_surface     # erg/s/cm2/A
flux_density=flux_density/0.07
master=spe.sum(0)
master=master/master[(wave>6583.)&(wave<6585.)].mean()*flux_density

waveHa,fluxHa=np.loadtxt('Halpha_model/Ha_60_14.dat',unpack=True)
waveHa=waveHa[::-1]
fluxHa=fluxHa[::-1]
fluxorig=fluxHa
# Planet radius = 2 R_J, emitting surface = 1% of planet surface
fluxHa=fluxHa*(2.*6.9911e7)**2*0.01/(113.4*3.0857e16)**2
# Normalization to get a line flux of 8.1e-16 erg/s/cm2 as for PDS 70b in Hashimoto et al. (2020)
dl = np.mean(np.diff(waveHa))
fluxHa = fluxHa / fluxHa.sum() / dl * 8.1e-16
#fluxHa=fluxHa/fluxHa.sum()/0.05123*8.1e-16
model=np.interp(wave,waveHa,fluxHa,left=0.,right=0.)

# SNR reference: SNR=45 per pixel of 0.82 km/s in 1800s for HARPS
# RISTRETTO: SNR=92 per pixel of 0.82 km/s in 3600s assuming total throughput of 3.0%
# To be improved: rebin master and model spectra to 0.82 km/s for better modeling of RON

contrast=500.
RON=6.
Fs=master/master[(wave>6583.)&(wave<6585.)].mean()*92**2*0.01/6562.82*3e5/0.82/contrast
Fp=model/master[(wave>6583.)&(wave<6585.)].mean()*92**2*0.01/6562.82*3e5/0.82
noise=np.random.standard_normal(len(Fs))*np.sqrt(Fs+Fp+RON**2)
simul=Fp+noise

plt.figure(1,figsize=(12.,7.))
plt.clf()
plt.plot(wave,master+model,'tab:orange')
plt.plot(wave,master,'tab:blue')
plt.xlim(6553.,6572.)
plt.ylim(2.7e-14,8.7e-14)
plt.xlabel('Wavelength (A)')
plt.ylabel('Flux Density (erg/s/cm2/A)')
plt.title('PDS 70 H-alpha Simulation')
plt.show()

plt.figure(2,figsize=(12.,7.))
plt.clf()
plt.plot(wave,master+model,'tab:orange')
plt.plot(wave,master,'tab:blue')
plt.xlim(6560.,6564.)
plt.ylim(4.5e-14,8.2e-14)
plt.xlabel('Wavelength (A)')
plt.ylabel('Flux Density (erg/s/cm2/A)')
plt.title('PDS 70 H-alpha Simulation')
plt.show()

plt.figure(3,figsize=(12.,7.))
plt.clf()
plt.plot(rv,simul,'tab:blue')
plt.xlim(-240.,90.)
plt.ylim(-25.,85.)
plt.xlabel('Radial Velocity (km/s)')
plt.ylabel('Flux (e-)')
plt.title('VLT-RISTRETTO Simulation of a Planet at 3-5 AU around PDS 70 (Texp=1 hour)')
plt.show()

plt.figure(4,figsize=(12.,7.))
plt.clf()
plt.plot(waveHa,FluxHa,'tab:blue')
plt.plot(waveHa,Fluxorig,'tab:blue')
plt.yscale('log')
#plt.xlim(-240.,90.)
#plt.ylim(-25.,85.)
plt.xlabel('Radial Velocity (km/s)')
plt.ylabel('Flux (e-)')
plt.title('VLT-RISTRETTO Simulation of a Planet at 3-5 AU around PDS 70 (Texp=1 hour)')
plt.show()

# Create output directory if it doesn't exist
outdir = "simulation_lovis"
os.makedirs(outdir, exist_ok=True)

# --- Save wave, master, and master+model ---
csv1 = os.path.join(outdir, "simul_lovis.csv")
np.savetxt(
    csv1,
    np.column_stack((wave, master, master + model)),
    delimiter=",",
    header="wavelength_A,master_flux,master_plus_model_flux",
    comments="",
)
print(f"Saved {csv1}")

# --- Save rv and simulated counts ---
csv2 = os.path.join(outdir, "simul_lovis_3-5au.csv")
np.savetxt(
    csv2,
    np.column_stack((rv, simul)),
    delimiter=",",
    header="rv_kms,simulated_flux_electrons",
    comments="",
)
print(f"Saved {csv2}")

