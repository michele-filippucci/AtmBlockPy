import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

import lib.BlockPlots as BP

climat1 = "/home/guest/work/michele/data/ERA5/processed/"+\
          "las_TM90_clim_djfm_northem_1979-2008.nc"
climat2 = "/home/guest/work/michele/data/ERA5/processed/"+\
          "lab_TM90_clim_djfm_northem_1979-2008.nc"

plot =  "/home/guest/work/michele/prog/plots/ensamble/t_student.png"

ds1 = xr.load_dataset(climat1)
ds2 = xr.load_dataset(climat2)
var1 = ds1["freq_variance"].values
var2 = ds2["freq_variance"].values
mean1 = ds1["freq_ensamble"].values
mean2 = ds2["freq_ensamble"].values
#print(var1,var2,mean1,mean2)

var1 = xr.where(var1==0,0.00001,var1)
var2 = xr.where(var2==0,0.00001,var2)

t_student = np.absolute(mean1-mean2)/(0.5*(var1**0.5 + var2**0.5))
ds1["freq_ERA5"][:,:,:] = t_student
print(np.amin(t_student),np.amax(t_student))

BP.PlotFreqZ500(ds=ds1,output=plot,mapcrs = ccrs.NorthPolarStereo(),freq="freq_ERA5",plot_title="T _ Student distribution")

sm = sum(sum(sum(np.absolute(t_student))))
shp = np.shape(t_student)
print(sm/(shp[0]*shp[1]*shp[2]))


