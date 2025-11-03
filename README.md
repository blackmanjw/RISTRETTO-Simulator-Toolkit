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
2. Downsample to spectral resolution of 140,000 (input is in the millions). <br><code>python resample.py PDS70c_BT-Settl-CIFIST-1300K-4logg.txt --R 140000</code><br> The --R parameter (desired resolution for resampling) is optional. The outputs have the extention 140000_orig.txt.
3. You can add the stellar H-alpha signal from PDS70 from HARPS with <code>python harps_pds70.py</code> to the the R~140000 resampled stellar spectra. This will create an output in the input directory that scales the stellar H-alpha signal from the HARPS data and combines it the BT-Settl to make PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.csv and a number of png plots. It also creates a comparison_Ha_blended.png plot.
4. Add the simulated H-alpha from Yuhiko Aoyama's model to the planet spectra. <br><code>python add-halpha.py 2MJ1612b_BT-Settl-CIFIST-1200K-3.5logg_140000.txt Ha_60_12.dat</code><br>You can add all four models in Halpha_models to the three planets (PDS70b+c, WISPIT2 and 2MJ1612b) by running <code>./addalpha.sh</code>
5. Convert to input for Pyechelle, eg. <br><code> python makecsv.py PDS70b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12.txt 100</code>. The 100 is optional and denotes the preferred separation for the off-axis case. If it is not included it will default to the on-axis case (coupling 0.45 for planet. For the star it will output versions both for a coupling of 0.45 and 3.5e-4). It is only required to use this for the star as the planet coupling is always the same (0.45). If set the code with draw on <code>ao-coupling/fibre_offaxis_vs_distance.txt</code> and interpolate to get the preferred value.
6. You can plot and overlay all the inputs to the above procedure (<code>python plotinput.py</code>) and also the outputs which will be the subsequent Pyecehlle inputs (the resulting .csv files) using the code <code>python plotoutput.py</code>
7. You can batch run all of the different separations at once using <br><code>python makescv_batch.py PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt</code> which will run makecsv for separations of 100-700mas at intervals of 50mas. This only needs to be run for the host stars. Below are the complete lines you need to run for the three objects (PDS70/WISPIT2/PDS70)<br>
<code>python makecsv_batch.py 2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000.txt<br>
python makecsv_batch.py PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt<br>
python makecsv_batch.py WISPIT2_star_BT-Settl-CIFIST-4400K-4logg_140000.txt<br>
python makecsv.py 2MJ1612b_BT-Settl-CIFIST-1200K-3.5logg_140000_Ha_60_12.txt<br>
python makecsv.py PDS70b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12.txt<br>
python makecsv.py PDS70c_BT-Settl-CIFIST-1300K-4logg_140000_Ha_60_12.txt<br>
python makecsv.py WISPIT2b_BT-Settl-CIFIST-1400K-4logg_140000_Ha_60_12.txt<br>
python makecsv.py 2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000.txt centralspaxel<br>
python makecsv.py PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt centralspaxel<br>
python makecsv.py WISPIT2_star_BT-Settl-CIFIST-4400K-4logg_140000.txt centralspaxel</code>
<br> You can also just skip all of this and run ./makeall.sh which will do everything (all planets and stars at different separations. You may need to make it executable <code>chmod +x makeall.sh</code>
8. Test the Pyechelle input with <br><code>python testinput.py PDS70_star_BT-Settl-CIFIST-4200K-5logg_140000.txt</code>. This adds noise and scales it to a set exposure time (by defauly 3600s). You can also run testinput_phoenix.py to run the same code by for a phoenix theoretical spectra. By default it is setup for PDS70.
9. Copy all the outputs to the Ubelix server so you can run Pyechelle. <br><code>scp output/*.csv jb23l046@submit03.unibe.ch:/storage/homefs/jb23l046/Simu_run/data/in/</code>
10. From Ubelix run <code>/storage/homefs/jb23l046/Simu_run/sbatch offaxis.sh</code> and <code>/storage/homefs/jb23l046/Simu_run/sbatch onaxis.sh</code>. This will run the corresponding onaxis.py and offaxis.py files which will batch process all the input files through Pyecehlle., In the on axis-case, we combine the star (0.45 coupling) in Fiber 1 with the Star (coupling 3.5e-4) and the planet (0.45 coupling) in Fiber 2. In the off-axis case we run the star through Fibre 1 (0.45 coupling) and separately run the star (with coupling determined by the off axis separation between 100-600mas) summed with the planet (0.45 coupling). This outputs to /data/out/offaxis and /data/out/onaxis.
11. Copy files from Ubelix to my PC. <code> scp -r jb23l046@submit03.unibe.ch:/storage/homefs/jb23l046/Simu_run/data/out/offaxis "/Users/odin/toolkit/pychelle_output"</code><br><code> scp -r jb23l046@submit03.unibe.ch:/storage/homefs/jb23l046/Simu_run/data/out/onaxis "/Users/odin/toolkit/pychelle_output"</code>.<br> To copy files (eg. the combined star and planet data/plot) from Ubelix, run <code>scp -r jb23l046@submit03.unibe.ch:/storage/homefs/jb23l046/Simu_run/data/in /Users/odin/toolkit/inputfromubelix/</code> You can run all of these at once with <code>./fromubelix.sh</code>
12. cp WISPIT2_star_0.45_3600s.fits, PDS70_star_0.45_3600s.fits, 2MJ1612_star_0.45_3600s.fits from offaxis to onaxis folders (this is for the star in the centre spaxel).
13. To run the analysis, run <code>python analysis.py onaxis/PDS70b_Ha_60_12_3600s</code>. This will generate some csv files and png files in the analysis subdirectory.
14. Run the interpretation code like <code>python interp2.py onaxis/PDS70b_star_plus_planet_Ha_60_12_3600s_fiber2_order117.csv --noise 6530 6550 6580 6605</code>. This will determine the signal-to-noise of the peak.

Command	Description
python analysis.py onaxis/PDS70b_Ha_80_14_3600s 2	Normal single-file spectrum
python analysis.py onaxis/PDS70b_Ha_80_14_3600s 2 all	Plot all 4 Ha variants for that object
python analysis.py plotall	Make 3 “all” plots (for WISPIT2b, PDS70b, 2MJ1612b) and then combine them into one figure with 3 subplots


16. You can run Lovis original approximate simulation with <code>python simulation_lovis.py</code>. This requires the data in the harps_pds70 and Halpha_model subdirectories.
#### Also

1. <code> python ao-coupling.py </code> is used to generate the adaptive optics coupling map for the different fibres (1,2, 3-7 averaged). For the off axis case we only consider Fibre 2. The output of this code is <code>ao-coupling/fibre_offaxis_vs_distance.txt</code>. This is read-in and interpolated when running makecsv.py with a specified separation.
2. Unit propagation of flux in makecsv.py
Input spectrum	erg / cm² / s / Å
Convert to meters	erg / cm² / s / m
Convert erg → J	J / cm² / s / m
Divide by photon energy	photons / cm² / s / m
cm² → m²	photons / m² / s / m
Multiply by star area	photons / s / m
Divide by distance²	photons / m² / s / m (at telescope)
Multiply by telescope area	photons / s / m
Multiply by efficiency	photons / s / m
Interpolate	photons / s / m
