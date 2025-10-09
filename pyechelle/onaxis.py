from pyechelle.simulator import Simulator
from pyechelle.spectrograph import ZEMAX
from pyechelle.sources import CSVSource
from pyechelle.telescope import Telescope
import numpy as np
from datetime import datetime, timedelta

base_path = '/storage/homefs/jb23l046/Simu_run/data/'

exp_time= 3600  # in seconds


for position in range(1):

    sim          = Simulator(ZEMAX("RISTRETTOtest18"))
    sim.set_ccd(1)
    #sim.set_cuda(True)
    sim.set_fibers([1,2]) # 1,2,3,4,5,6,7
    sim.set_sources([
        CSVSource(base_path + f"in/2MJ1612_star_BT-Settl-CIFIST-3900K-4logg_140000_fiber1_0.45_central.csv", list_like=False,wavelength_units='nm',flux_units='ph/s/AA'),
        CSVSource(base_path + f"in/fiber_2.csv", list_like=False,wavelength_units='nm',flux_units='ph/s/AA')])
    #sim.set_atmospheres([True], sky_calc_kwargs={'airmass': airmass[position],'observatory':'paranal'})
    sim.set_exposure_time(exp_time)
    sim.set_bias(250)
    sim.set_read_noise(3)
    sim.set_output(base_path + 'out/onaxis/'+f'test0810_{exp_time}_s.fits', overwrite=True)
    sim.run()
