import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

np.set_printoptions(precision=2,threshold=np.inf)

fn = "/home/guest/work/michele/data/ERA5/ERA5_pIB_daily_djfm_northem_1979-2019.nc"
fn2 = "/home/guest/work/michele/data/ERA5/ERA5_precipitation_day_djfm_northem_r144x73_1979-2019_bis.nc"
dat_nam1 = "zg"
dat_nam2 = "pIB_boolean"
dat_nam3 = "tp"
#day = "1979-01-01T09:00:00"
plot_title = 'Z500 composites on European Blocking vs precipitation anomaly'
longitude =10.
latitude = 55.


# Import Data
ds = xr.open_dataset(fn)
ds2 = xr.open_dataset(fn2)
lats = ds.lat.data
lons = ds.lon.data
dat1_zero = ds[dat_nam1].loc[:,5e+04,:,:]/9.80665
dat2_zero = ds[dat_nam2].loc[:,5e+04,:,:]
dat3 = ds2[dat_nam3].loc[:,:,:]

counter = 0
meanT = dat3.values.mean(0)

for time in ds.time.values:
    if dat2_zero.loc[time,latitude,longitude].values == 0:
        dat1_zero.loc[time,:,:] = 0
        dat3.loc[time,:,:] = 0
        counter += 1

dat1_zero = dat1_zero.values
dat3 = dat3.values
#dat1 = xr.where(x == 0, dat1_zero, 0) #set to zero the array1 when array2 == 0
#dat2 = dat2_zero.mean(0)*100#using time, not zero because the type is xarray
dat1 = np.ma.masked_equal(dat1_zero,0)#mask array1 when ==0
dat3 = np.ma.masked_equal(dat3,0)
print(dat1.shape)
dat1 = dat1.mean(0)# - dat1_zero.mean(0)
dat3 = dat3.mean(0) - meanT #temp anomaly
#dat1 = dat1_zero.values.mean(0) #uncomment if you want mean over ERA5 period


#check output wanted
#print(x.shape)
perc = (1-(counter/dat2_zero.shape[0]))*100
print(perc)
print(dat1_zero.shape)

# Set up the projection that will be used for plotting the map
mapcrs = ccrs.NorthPolarStereo()
# Set up the projection of the data; if lat/lon then PlateCarree is what you want
datacrs = ccrs.PlateCarree()

#plot map subplot
ax = plt.subplot(111, projection=mapcrs)
ax.set_extent([-180, 175.5, 0, 90], ccrs.PlateCarree())
ax.coastlines()

#plot contour
cs_range = np.arange(5000,6000,50)
#cb_range = np.arange(0,1,0.05)

cs = ax.contourf(lons, lats,dat3,cmap="coolwarm", transform=datacrs) #first ploot (float)
cs_b = ax.contour(lons,lats,dat1, cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#second plot for contours
plt.clabel(cs_b, colors  ="black", fontsize ="smaller", inline ="false")
#cb = ax.contour(lons,lats, dat2,cmap = "Greys", transform=datacrs)#blocking plot

# Make some nice titles for the plot (one right, one left)
plt.title(plot_title, loc='left')
#plt.title("% events over the period: " + str(perc)[0:5], loc="right")
plt.plot(longitude,latitude,"or",transform=datacrs)
plt.colorbar(cs, orientation='horizontal', pad=0, label = "Temperature (K)" , aspect=40)
#plt.colorbar(cb, orientation='horizontal', label = 'blocking freq (%)', aspect = 50)
#plt.colorbar.label("% over the period")
# Show image
plt.show()
