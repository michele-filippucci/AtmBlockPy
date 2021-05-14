import numpy as np
import xarray as xr

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
fn2 = "/home/guest/work/michele/data/ERA5/T500/ERA5_T500_day_djfm_northem_r144x73_1979-2019.nc"

fn_out = "/home/guest/work/michele/data/ERA5/processed/ERA5_pIB_daily_djfm_northem_1979-2019.nc"

img_out = "/home/guest/work/michele/prog/plots/zg500composite_tempanomaly.png"
"""
#try BlockTools
pIB = BlockTools()
pIB.read(fn)
print(pIB.boolean_pIB(fn_out,freq_also = True))
"""
plot = BlockPlots("ZG500 composite vs temp anomaly")
plot.read_main(fn_out)
plot.read_additional_ds(fn2)
#try composite
print(plot.PlotCompositeZ500(output = img_out,additional = "ta",point_coords=[55.0,10.0]))
