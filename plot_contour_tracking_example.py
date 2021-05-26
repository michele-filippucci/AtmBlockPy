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

fn = "/home/guest/work/michele/data/ERA5/processed/"+\
         "ERA5_pIB_daily_djfm_northem_1995-2000.nc"
fn_out = "/home/guest/work/michele/data/ERA5/processed/"+\
         "ERA5_pIB_tracked_daily_djfm_northem_1995-2000.nc"
img_out = "/home/guest/work/michele/prog/plots/"+\
          "contour_tracking/"+\
          "TrackingExample_1997-01-01_wpers4_50%b.png"
print("Starting job")
#try contour tracking
contourpIB = BlockTools()
contourpIB.read(fn)
print("Data correctly read")
contourpIB.ContourTracking2D(fn_out,pers=4)
print("Data created")
plot = BlockPlots("tracking 1997-01-01T09:00:00")
plot.read_main(fn_out)
print(plot.PlotTracking(output = img_out,\
      starting_day="1997-01-01T09:00:00.000000000"))
print("Plot produced")

print("--- %s seconds ---" % (time.time() - start_time))

