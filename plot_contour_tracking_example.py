import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

#my own classes
import sys
sys.path.append("/home/guest/work/michele/prog/GitPROG/AtmBlockPy/lib")
from BlockTools import BlockTools
import BlockPlots as BP
import BlockTools as BT

#list = ["01-01"] #,"01-15","02-01","02-15","03-01","03-15"]
import time
start_time = time.time()

fn = "/home/guest/work/michele/data/ERA5/Z500/"+\
         "ERA5_daily_jja_northem_2003.nc"
fn_out = "/home/guest/work/michele/data/ERA5/processed/"+\
         "ERA5_pIB_tracked_daily_jja_northem_2003.nc"

#try contour tracking

pIB = BlockTools()
pIB.read(fn)
print("Data correctly read")
ds = pIB.TM90(data_return=True)
contourpIB = BlockTools()
contourpIB.load_data(ds)
contourpIB.ContourTracking2D(fn_out,pers=0)
print("Data created")

img_out = "/home/guest/work/michele/prog/plots/"+\
        "CNN/"+\
        "CaseStudy.png"
print("Starting job")
ds = xr.load_dataset(fn_out)
print(BP.PlotTracking(ds=ds,output = img_out,\
      starting_day="2003-"+"07-28"+"T09:00:00.000000000",extent=[-10,40,30,75]))
print("Plot produced")
  
print("--- %s seconds ---" % (time.time() - start_time))

	
