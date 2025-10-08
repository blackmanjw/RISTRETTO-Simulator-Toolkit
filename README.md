# toolkit

## baraffe_isochrone.py
Output interpolated values for T_eff and logg from an input age and mass using the COND03_models.txt
This file is the data discussed in: Baraffe, Chabrier, Barman, Allard, Hauschildt, 2003, A&A, accepted "Evolutionary models for cool brown dwarfs and extrasolar giant planets. The case of HD 209458"
OR https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/ BHAC15_iso.2mass which are newer isochrones that go to higher masse, Models from Baraffe, Homeier, Allard, Chabrier 2015, A&A,577, 42
"New evolutionary models for pre-main sequence and main sequence low-mass stars down to the hydrogen-burning limit"

NOTE: Isochrones in different filters will be regularly added on this site.

The file BHAC15_tracks+structure include tracks and some information on the internal structure
(mass of the radiative core, gyration radii, central temperature and densities, etc...).

## RISTRETTO simulation procedure

1. Download desired BT-Settl model (set temperature and log(g)) from https://svo2.cab.inta-csic.es/theory/newov2/index.php?models=bt-settl-cifist. If necessary use the baraffe_isochrone tool to determine the Teff and log(g) from the mass and age.
2. Downsample to spectral resolution of 140,000 (input is in the millions). <br><code>python resample.py PDS70c_BT-Settl-CIFIST-1300K-4logg.txt --R 140000</code><br> The --R parameter (desired resolution for resampling) is optional.
3. Add the simulated H-alpha from Yuhiko Aoyama's model. <br><code>python add-halpha.py 2MJ1612b_BT-Settl-CIFIST-1200K-3.5logg_140000.txt Ha_60_12.dat</code>
4. Convert to input for Pyechelle, eg. <br><code> python makecsv.py PDS70b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12.txt 100</code>. The 100 is optional and denotes the preferred separation for the off-axis case. If it is not included it will default to the on-axis case (coupling 0.45 for planet, 3.5e-4 for star). It is only required to use this for the star as the planet coupling is always the same (0.45). If set the code with draw on <code>ao-coupling/fibre_offaxis_vs_distance.txt</code> and interpolate to get the preferred value. 
5. You can plot and overlay all the inputs to the above procedure (<code>python plotinput.py</code>) and also the outputs which will be the subsequent Pyecehlle inputs (the resulting .csv files) using the code <code>python plotoutput.py</code>
6. You can batch run all of the different separations at once using <br><code>python makescv_batch.py PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt</code> which will run makecsv for separations of 100-700mas at intervals of 50mas. This only needs to be run for the host stars. Below are the complete lines you need to run for the three objects (PDS70/WISPIT2/PDS70)<br>
<code>python makecsv_batch.py 2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000.txt<br>
python makecsv_batch.py PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt<br>
python makecsv_batch.py WISPIT2_star_BT-Settl-CIFIST-4400K-4logg_140000.txt<br>
python makecsv.py 2MJ1612b_BT-Settl-CIFIST-1200K-3.5logg_140000_Ha_60_12.txt<br>
python makecsv.py PDS70b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12.txt<br>
python makecsv.py PDS70c_BT-Settl-CIFIST-1300K-4logg_140000_Ha_60_12.txt<br>
python makecsv.py WISPIT2b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12.txt</code>
6. Test the Pyechelle input with <br><code>python testinput.py WISPIT2b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12_fiber1.csv 1800</code>. This adds noise and scales it to a set exposure time (here 1800s). Adding exposure time is optional. The default is 3600s (1 hour).
7. Copy all the outputs to the Ubelix server so you can run Pyechelle. <br><code>scp output/*.csv jb23l046@submit03.unibe.ch:/storage/homefs/jb23l046/Simu_run/data/in/</code>

#### Also

1. <code> python ao-coupling.py </code> is used to generate the adaptive optics coupling map for the different fibres (1,2, 3-7 averaged). For the off axis case we only consider Fibre 2. The output of this code is <code>ao-coupling/fibre_offaxis_vs_distance.txt</code>. This is read-in and interpolated when running makecsv.py with a specified separation.
