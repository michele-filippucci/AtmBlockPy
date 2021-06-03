import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

#my own class
from lib.BlockTools import BlockTools
from lib.BlockPlots import BlockPlots
import lib.BlockTools as BT

from cartopy.examples.waves import sample_data

import time
start_time = time.time()

fn = "/home/guest/work/michele/data/ERA5/Z500/"+\
     "ERA5_Z500_day_djfm_r144x73_500hPa_northem_1979-2019.nc"
fn_out = "/home/guest/work/michele/data/ERA5/processed/"+\
         "ERA5_pIB_daily_wfilters_djfm_northem_1979-2019.nc"
img_out = "/home/guest/work/michele/prog/plots/filters/"+\
          "Freq1979-2019wlong_filter.png"


#BlockTools
pIB = BlockTools()
pIB.read(fn)
print(pIB.boolean_pIB(fn_out,freq_also = True,mer_gradient_filter = False))

#BlockPlots
freqPlot = BlockPlots("1979-2019 Blocking Freq w long_filter")
freqPlot.read_main(fn_out)
print(freqPlot.PlotFreqZ500(img_out))

print("--- %s seconds ---" % (time.time() - start_time))
