import numpy as np
import xarray as xr

class BlockTools(object):
    num_of_BlockTools = 0

    def __init__(self):
        BlockTools.num_of_BlockTools += 1

    def read(self,filename):
        self.dataset = xr.load_dataset(filename)

    def boolean_pIB(self,fn_out = "",data_return = False,freq_also = False):
        if fn_out=="" and data_return==False:
            string = "Specify the kind of output you want"
            print(string)
            return 0
        #checking if dataset is right
        try:
            zg = self.dataset["zg"]
            string = "data successfully received"
        except:
            string = "zg variable was not found.\n\ Hint: use read() to load data."
            print(string)
            return 0

        #____CHECK GEOP HEIGHT____
        #ERA5 dataset uses geopotential, not height
        if zg.values[0,0,0,0] > 10000:
        	zg = zg/9.80665

        #.values gives tuples
        #compute GHGS/GHGN
        GHGS = (+ zg.loc[:,:,30:75,:].values - zg.loc[:,:,15:60,:].values)/15.0
        GHGN = (- zg.loc[:,:,30:75,:].values + zg.loc[:,:,45:90,:].values)/15.0

        #look for grid points where the conditions for pIB are satisfied
        #using where function from xarray
        #first term is GHGN condition, second term is GHGS condition. Tuples are multiplied
        #(no matrix product)
        TuplepIB = xr.where(GHGN < -10.0, 1.0, 0.0) * xr.where(GHGS > 0., 1.0 , 0.0)

        #define XArray using the costructor for an appropriate output
        times = zg.coords["time"].values
        plev = zg.coords["plev"].values
        lon = zg.coords["lon"].values
        lat = zg.coords["lat"].values
        pIB = xr.DataArray(0,coords=[times,plev,lat,lon],dims = zg.dims)
        pIB.loc[:,:,30:75,:] = TuplepIB
        self.dataset = self.dataset.assign(pIB_boolean = pIB)
        if freq_also == True:
            pIB_f = sum(pIB)*100/pIB.values.shape[0]
            self.dataset = self.dataset.assign(pIB_frequencies = pIB_f)
        if data_return == False:
            self.dataset.to_netcdf(fn_out)
        if data_return == True:
            print(string)
            return self.dataset
        else:
            return 0
