import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

#my own class
import lib.BlockTools as BT
import lib.BlockPlots as BP
from cartopy.examples.waves import sample_data

import time
start_time = time.time()

fn = "/home/guest/work/michele/data/ERA5/processed/ERA5_pIB_daily_djfm_northem_1979-2019.nc"
#fn2 = "/home/guest/work/michele/data/ERA5/T500/"+\
#       "ERA5_T500_day_djfm_northem_r144x73_1979-2019.nc"
fn2 = "/home/guest/work/michele/data/ERA5/precipitation/"+\
      "ERA5_pr_day_r144x73_djfm_northem_1979-2019.nc"
img_out = "/home/guest/work/michele/prog/plots/composites/"+\
          "ZG500CompositePrecipitation55N10E_1979-2019_NorthPolarStereo.png"
print("Starting job")
#try composite precipitation
ds = xr.load_dataset(fn)
add_ds = xr.load_dataset(fn2)
print("Data correctly read")
print(BP.PlotCompositeZ500(ds = ds,
                             add_ds = add_ds,
                             output = img_out,\
                             additional ="pr",\
                             mapcrs = ccrs.NorthPolarStereo(),\
                             point_coords=[55.0,10.0],\
                             colorbar_label = "Precipitation anomaly (mm)",\
                             extent=[-180,180,0,90],
                             plot_title="ZG500 composite vs prec anomaly, 55 N, 10 E"))
print("Plot produced")
print("--- %s seconds ---" % (time.time() - start_time))
