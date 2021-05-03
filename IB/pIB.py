"""

Punctual Instantaneus Blocking Detection

This program performs the detection of the instantaneus blocking in each
grid point of the dataset examined. It outputs a new dataset containing two new
variables:

- pIB_boolean = big matrix containing 1,0 values. 1 for blocking event, 0 for
  not blocking event
- pIB_frequencies = matrix without time dimension indicating for what percentual of
  time a grid point has been marked with 1 in the matrix above.

"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

fn = "/home/guest/work/michele/data/ERA5/ERA5_Z500_day_r144x73_500hPa_northem_1979-2019.nc" #nome del file di partenza
fn_out = "/home/guest/work/michele/data/ERA5/ERA5_pIB_daily_northem_1979-2019.nc"
#____CDO FILTERS____
#/

#____IMPORT DATA____
#store zg in a xarray
data = xr.load_dataset(fn)
zg = data["zg"]
print(data)

#____CHECK GEOP HEIGHT____
#ERA5 dataset uses geopotential, not height
if zg.values[0,0,0,0] > 10000:
	zg = zg/9.80665

#____DEF punctual IB____
def pIB(ZG):

	#.values gives tuples
	#compute GHGS/GHGN
	GHGS = (+ ZG.loc[:,:,30:75,:].values - ZG.loc[:,:,15:60,:].values)/15.0
	GHGN = (- ZG.loc[:,:,30:75,:].values + ZG.loc[:,:,45:90,:].values)/15.0

	#look for grid points where the conditions for pIB are satisfied
	#using where function from xarray
	#first term is GHGN condition, second term is GHGS condition. Tuples are multiplied
	#(no matrix product)
	TuplepIB = xr.where(GHGN < -10.0, 1.0, 0.0) * xr.where(GHGS > 0., 1.0 , 0.0)

	#define XArray using the costructor for an appropriate output
	times = ZG.coords["time"].values
	plev = ZG.coords["plev"].values
	lon = ZG.coords["lon"].values
	lat = ZG.coords["lat"].values
	pIB = xr.DataArray(0,coords=[times,plev,lat,lon],dims = ZG.dims)
	pIB.loc[:,:,30:75,:] = TuplepIB
	return pIB


#CHECK RESULTS
pIB_b = pIB(zg) #pIB boolean

#add punctual IB boolean variable to dataset
data = data.assign(pIB_boolean = pIB_b)

#add punctual IB frequencies variable to dataset
pIB_f = sum(pIB_b)*100/pIB_b.values.shape[0]
data = data.assign(pIB_frequencies = pIB_f)

print(data)
data.to_netcdf(fn_out)
