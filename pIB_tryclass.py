import numpy as np
import xarray as xr

#my own class
from BlockTools import BlockTools
from BlockPlots import BlockPlots

fn = "/home/guest/work/michele/data/ERA5/Z500/ERA5_Z500_day_djfm_r144x73_500hPa_northem_1979-2019.nc"
fn_out = "/home/guest/work/michele/data/ERA5/processed/ERA5_pIB_daily_djfm_northem_1979-2019.nc"

img_out = "/home/guest/work/michele/prog/local/plots/zg500fieldblockingevent.png"
"""
#try BlockTools
pIB = BlockTools()
pIB.read(fn)
print(pIB.boolean_pIB(fn_out,freq_also = True))
"""
plot = BlockPlots("ZG500 field - Blocking Event")
plot.read_main(fn_out)
print(plot.PlotEventZ500(output = img_out,day= "1979-01-01T09:00:00"))
