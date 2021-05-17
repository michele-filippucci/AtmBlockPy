import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

#my own class
from BlockTools import BlockTools
from BlockPlots import BlockPlots


from cartopy.examples.waves import sample_data

"""
x, y, z = sample_data((10, 20))

print(x)
print(y)
print(z)
"""

fn = "/home/guest/work/michele/data/ERA5/Z500/ERA5_Z500_day_djfm_r144x73_500hPa_northem_1979-2019.nc"
#fn2 = "/home/guest/work/michele/data/ERA5/T500/ERA5_T500_day_djfm_northem_r144x73_1979-2019.nc"
fn2 = "/home/guest/work/michele/data/ERA5/precipitation/ERA5_precipitation_day_djfm_northem_r144x73_1979-2019.nc"

fn_out = "/home/guest/work/michele/data/ERA5/processed/ERA5_pIB_daily_djfm_northem_1979-2019.nc"
fn_out2 = "/home/guest/work/michele/data/ERA5/processed/ERA5_pIB_tracked_daily_djfm_northem_1979-2019.nc"
img_out = "/home/guest/work/michele/prog/plots/Tracking.png"

#try contour tracking
contourpIB = BlockTools()
contourpIB.read(fn_out)
#ds = contourpIB.ContourTracking(data_return = True)
#print(ds["pIB_tracked"])
contourpIB.ContourTracking(fn_out2)

plot = BlockPlots("tracking 1980-01-01T09:00:00 - 1980-01-15T09:00")
plot.read_main(fn_out2)
#plot.read_additional_ds(fn2)
print(plot.PlotTracking(output = img_out,days=["1979-01-01T09:00:00.000000000","1979-01-16T09:00:00.000000000"]))
"""
#try BlockTools
pIB = BlockTools()
pIB.read(fn)
print(pIB.boolean_pIB(fn_out,freq_also = True))
"""

"""
#try composite temperature
plot = BlockPlots("ZG500 composite vs temp anomaly, 55 N, 10 E")
plot.read_main(fn_out)
plot.read_additional_ds(fn2)
#try composite
print(plot.PlotCompositeZ500(output = img_out,additional = "ta",point_coords=[55.0,10.0]))
"""
"""
#try composite precipitation
plot = BlockPlots("ZG500 composite vs prec anomaly, 55 N, 10 E")
plot.read_main(fn_out)
plot.read_additional_ds(fn2)
print(plot.PlotCompositeZ500(output = img_out, additional ="tp",mapcrs = ccrs.PlateCarree(),point_coords=[55.0,10.0],colorbar_label = "Precipitation (mm)",extent=[-70,70,20,70]))
"""

