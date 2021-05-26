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

fn = "/home/guest/work/michele/data/ERA5/processed/ERA5_pIB_daily_djfm_northem_1979-2019.nc"
fn2 = "/home/guest/work/michele/data/ERA5/T500/"+\
       "ERA5_T500_day_djfm_northem_r144x73_1979-2019.nc"
#fn2 = "/home/guest/work/michele/data/ERA5/precipitation/"+\
#      "ERA5_precipitation_day_djfm_northem_r144x73_1979-2018.nc"
img_out = "/home/guest/work/michele/prog/plots/composites/"+\
          "ZG500CompositeTemperature55N10E_1979-2019_NorthPolarStereo.png"

print("Starting job")
#try composite precipitation
plot = BlockPlots("ZG500 composite vs temp anomaly, 55 N, 10 E")
plot.read_main(fn)
plot.read_additional_ds(fn2)
print("Data correctly read")
print(plot.PlotCompositeZ500(output = img_out,\
                             additional ="ta",\
                             mapcrs = ccrs.NorthPolarStereo(),\
                             point_coords=[55.0,10.0],\
                             colorbar_label = "Temperature anomaly (K)",\
                             extent=[-180,180,0,90]))
print("Plot produced")
print("--- %s seconds ---" % (time.time() - start_time))
