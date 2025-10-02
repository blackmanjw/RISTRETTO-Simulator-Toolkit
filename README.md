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
4. Convert to input for Pyechelle, eg. <br><code> python makecsv.py PDS70b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12.txt</code>
5. You can plot and overlay all the inputs to the above procedure (<code>python plotinput.py</code>) and also the outputs which will be the subsequent Pyecehlle inputs (the resulting .csv files) using the code <code>python plotoutput.py</code>
6. Test the Pyechelle input with <br><code>python testinput.py WISPIT2b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12_fiber1.csv 1800</code>. This adds noise and scales it to a set exposure time (here 1800s). Adding exposure time is optional. The default is 3600s (1 hour).
7. Copy all the outputs to the Ubelix server so you can run Pyechelle. <br><code>scp output/*.csv jb23l046@submit03.unibe.ch:/storage/homefs/jb23l046/Simu_run/data/in/</code
