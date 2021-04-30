"""
Short program to confront outputs
-This is written for a specific dataset and is not versatile

"""

#piccola variazione

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

fn1 = "/home/guest/work/michele/data/ERA5/ERA5_pIB_daily_djfm_northem_2015-2019.nc" #mio file
fn2 = "/home/guest/work/michele/data/CLIM_CRSCHK/D12_Full_ERA5_2015_2019_DJFM.nc"
fn3 = "/home/guest/work/michele/data/CLIM_CRSCHK/D12_Clim_ERA5_2015_2019_DJFM.nc"

data1= xr.load_dataset(fn1,decode_times=False)
data2= xr.load_dataset(fn2,decode_times=False)
data3= xr.load_dataset(fn3,decode_times=False)

pIB1= data1.pIB_boolean.values
pIB2= data2.InstBlock.values

pIB1_f= data1.pIB_frequencies.values
pIB3_f= data3.InstBlock.values

pIBdiff = data1.pIB_boolean.values - data2.InstBlock.values
pIB_fdiff = pIB1_f - pIB3_f

print(sum(sum(sum(sum(data1.pIB_boolean.values)))))
print(sum(sum(sum(sum(data2.InstBlock.values)))))
print(sum(sum(sum(sum(pIBdiff)))))

print("\n\n")

print(sum(sum(sum(pIB1_f))))
print(sum(sum(sum(sum(pIB3_f)))))
print(sum(sum(sum(sum(pIB_fdiff)))))
