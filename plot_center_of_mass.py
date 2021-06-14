import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil
import matplotlib.pyplot as plt

#my own class
from lib.BlockTools import BlockTools
import lib.BlockPlots as BP
import lib.BlockTools as BT

fn = "/home/guest/work/michele/data/ERA5/processed/"+\
         "ERA5_pIB_daily_djfm_northem_1995-2000.nc"
fn_out = "/home/guest/work/michele/data/ERA5/processed/"+\
         "ERA5_pIB_tracked_daily_djfm_northem_1995-2000.nc"

print("____PLOT: 01 starting____")
img_out = "/home/guest/work/michele/prog/plots/"+\
          "contour_tracking/series/"
print("Starting job")
ds = xr.load_dataset(fn_out)
tuple = ds["pIB_tracked"].values
tuple = tuple[0:10,0,:,:]
tuple = BT.OrderIndexes(tuple)
x,y = BT.CenterofMass(tuple,2)
print(len(x),len(y))
plt.plot(x, y)
plt.show()
