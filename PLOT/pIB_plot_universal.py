import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

fn = "/home/guest/work/michele/data/ERA5/ERA5_pIB_daily_northem_1979-2019.nc"
dat_nam1 = "zg"
dat_nam2 = "pIB_boolean"
day = "1979-01-01T09:00:00"
plot_title = '500-hPa ERA5 zg ' + day


# Import Data
ds = xr.open_dataset(fn)
lats = ds.lat.data
lons = ds.lon.data
dat1_zero = ds[dat_nam1]
dat2_zero = ds[dat_nam2]
dat1 = dat1_zero.loc[day,5e+04,:,:].values/9.80665
dat2 = dat2_zero.loc[day,5e+04,:,:].values

# Set up the projection that will be used for plotting the map
mapcrs = ccrs.PlateCarree()
# Set up the projection of the data; if lat/lon then PlateCarree is what you want
datacrs = ccrs.PlateCarree()

#plot map subplot
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([-15, 70, 30, 80], ccrs.PlateCarree())
ax.coastlines()

#plot contour
cs_range = np.arange(4600,6100,40)
cb_range = range(0,1)

cs = ax.contourf(lons, lats,dat1, cs_range,cmap="jet", transform=datacrs) #first ploot (float)
cs_b = ax.contour(lons, lats,dat1, cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#second plot for contours
plt.clabel(cs_b, colors  ="black", fontsize ="smaller", inline ="false")
cb = ax.contour(lons,lats, dat2, cb_range,cmap = "flag", transform=datacrs)#blocking plot

# Make some nice titles for the plot (one right, one left)
plt.title(plot_title, loc='left')
plt.colorbar(cs, orientation='horizontal', pad=0, label = "geopotential height" , aspect=50)
#plt.colorbar.label("% over the period")
# Show image
plt.show()
