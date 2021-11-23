"""
This example shows how to visualize the output of tracking. A dataset with a "_tracked" attribute is imported
and a plot of a given event (identified by the label) is performed.
In this example I use a climate SPHINX dataset.
"""

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

#my own classes
import sys
sys.path.append("/home/guest/work/michele/prog/GitPROG/AtmBlockPy/lib")
import BlockTools as BT
import BlockPlots as BP

#_______

label = 32
ens = 4 #ensemble member

#_______


ERA5 = "data/ERA5/processed/"+\"ERA5_pIB_tracked_daily_djfm_northem_1979-2019.nc"
lab = "SPHINX/processed/DAY_lab"+str(ens)+"_TM90_Tracked_clim_djfm_northem_1979-2008.nc"
las = "SPHINX/processed/DAY_las"+str(ens)+"_TM90_Tracked_clim_djfm_northem_1979-2008.nc"
mab = "SPHINX/processed/DAY_mab"+str(ens)+"_TM90_Tracked_clim_djfm_northem_1979-2008.nc"
mas = "SPHINX/processed/DAY_mas"+str(ens)+"_TM90_Tracked_clim_djfm_northem_1979-2008.nc"

#_______

dataset = mab

#_______

img_out = "/home/guest/work/michele/prog/plots/events/event.png"

print("____PLOT label " + str(label) + " starting____")
print("Starting job")
ds = xr.load_dataset(dataset)
print(BP.PlotLabel(ds=ds,output = img_out,shuffle=False,label = label))
print("Plot produced")





