import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
#import metpy.calc as mpcalc
#from metpy.units import units
import numpy as np
#from scipy.ndimage import gaussian_filter
import xarray as xr

fn = "/home/guest/work/michele/data/ERA5/ERA5_pIB_daily_northem_1979-2019.nc"

# Import Data
ds = xr.open_dataset(fn)
lats = ds.lat.data
lons = ds.lon.data
freq_z = ds["pIB_frequencies"].values
freq = freq_z[0,:,:]

# Set up the projection that will be used for plotting the map
mapcrs = ccrs.PlateCarree()
# Set up the projection of the data; if lat/lon then PlateCarree is what you want
datacrs = ccrs.PlateCarree()

#plot map subplot
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([-180, +180, 0, 90], ccrs.PlateCarree())
ax.coastlines()

#plot contour
cs = ax.contourf(lons, lats, freq,cmap=plt.cm.BuPu,transform=datacrs)

# Make some nice titles for the plot (one right, one left)
plt.title('500-hPa ERA5 Inst Puncutal Blocking frequency 2015/2019 djfm', loc='left')
plt.colorbar(cs, orientation='horizontal', pad=0, label = "% over the period" , aspect=50)
#plt.colorbar.label("% over the period")
# Show image
plt.show()
