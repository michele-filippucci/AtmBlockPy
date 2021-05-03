import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

fn = "/home/guest/work/michele/data/ERA5/ERA5_pIB_daily_northem_1979-2019.nc"
dat_nam1 = "zg"
dat_nam2 = "pIB_boolean"
day = "1979-01-01T09:00:00"
plot_title = 'Mean ZG during blocking events vs blocking frequency'


# Import Data
ds = xr.open_dataset(fn)
lats = ds.lat.data
lons = ds.lon.data
dat1_zero = ds[dat_nam1].loc[:,5e+04,:,:]/9.80665
dat2_zero = ds[dat_nam2].loc[:,5e+04,:,:]


dat1 = xr.where(dat2_zero == 1, dat1_zero, 0) #set to zero the array1 when array2 == 0
dat2 = dat2_zero.values.mean(0)*100#using time, not zero because the type is xarray
#dat1 = xr.where(dat2 > 5, dat1, 0)
dat1 = np.ma.masked_equal(dat1,0)#mask array1 when ==0
#dat1 = dat1[0,:,:]
dat1 = dat1.mean(0)
#check output wanted
print(dat1.shape)
print(dat1_zero.shape)

# Set up the projection that will be used for plotting the map
mapcrs = ccrs.PlateCarree()
# Set up the projection of the data; if lat/lon then PlateCarree is what you want
datacrs = ccrs.PlateCarree()

#plot map subplot
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([-180, +180, 30, 75], ccrs.PlateCarree())
ax.coastlines()

#plot contour
cs_range = np.arange(5000,6000,20)
#cb_range = np.arange(0,1,0.05)

cs = ax.contourf(lons, lats,dat1, cs_range,cmap="jet", transform=datacrs) #first ploot (float)
cs_b = ax.contour(lons, lats,dat1, cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#second plot for contours
plt.clabel(cs_b, colors  ="black", fontsize ="smaller", inline ="false")
cb = ax.contour(lons,lats, dat2,cmap = "Blues", transform=datacrs)#blocking plot

# Make some nice titles for the plot (one right, one left)
plt.title(plot_title, loc='left')
plt.colorbar(cs, orientation='horizontal', pad=0, label = "geopotential height (m)" , aspect=50)
plt.colorbar(cb, orientation='horizontal', label = 'blocking freq (%)', aspect = 50)
#plt.colorbar.label("% over the period")
# Show image
plt.show()
