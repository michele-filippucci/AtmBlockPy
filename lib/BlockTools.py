import numpy as np
import xarray as xr
from scipy.ndimage.measurements import label
#from skimage.morphology import diamond

np.set_printoptions(precision=2,threshold=np.inf)

def GetIndex(ds,coord="",key=""):
  index = 0
  for x in ds[coord].values:
    if str(x) == key:
      break
    index += 1
  return index

def OrderIndexes(arr):
  boolarr = arr > 0
  newarr = arr[boolarr]
  newarr = np.unique(np.sort(newarr))
  newval = range(1,len(newarr)+1)
  for i in range(0,len(newarr)):
    boolarr = arr == newarr[i]
    arr = xr.where(boolarr, newval[i],arr)
  return arr

class BlockTools(object):
  num_of_BlockTools = 0

  def __init__(self):
    BlockTools.num_of_BlockTools += 1

  def read(self,filename):
    self.dataset = xr.load_dataset(filename)


  """
  This function is capable of creating an nc file identical (located in fn_out)
  to the BlockTools.dataset with additional attributes:
  pIB_boolean and (when freq_also == True) pIB_frequencies
  As an alternative you can change the flag "data_return" and the function
  will return a dataset object from the class xarray containing the same
  additional attributes
  """

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

  """
  Contour Tracking
  This function takes a .nc file containing the pIB_boolean attribute from the function
  boolean_pIB and creates a new .nc file containing an additional attribute called 
  pIB_tracked which is zero when blocking is not occuring and is n when the nth blocking
  event is occuring.
  """
  def ContourTracking(self,fn_out = "",var_name = "pIB_boolean",data_return = False):
    if fn_out=="" and data_return==False:
      string = "Specify the kind of output you want"
      print(string)
      return 0
    try:
      pIB_boolean = self.dataset[var_name]
    except:
      print("Error Code 1: dataset not valid. The variable " + var_name + " cannot be found")
      return 1
    arr = pIB_boolean.values[:,0,:,:] 

    #filtra circa 9 punti griglia

    #label method from scipy.ndimage.measurements is used
    #structure = np.ones((3,3,3))
    structure = [[[0,0,0],[0,1,0],[0,0,0]],\
                 [[0,1,0],[1,1,1],[0,1,0]],\
                 [[0,0,0],[0,1,0],[0,0,0]]] #this matrix defines what is defined as neighbour
    #neighbour points are labeled with the same sequential number
    arr,ncomponents=label(arr,structure)
    #print(ncomponents)
    #applying some filters
    for t in np.arange(0,len(self.dataset.time.values-1)):
      bool = arr[t,:,:] > 0
      list = arr[t,bool]
      #print(np.unique(list))
      for l in np.unique(list):
        boolarr = arr[t,:,:] == l
        n = np.count_nonzero(boolarr)
        #print(n)
        #filtering cluster dimension
        if n < 9:
          #print("filtered")
          arr[t,:,:] = xr.where(boolarr, 0,arr[t,:,:])
          #print(np.count_nonzero(arr[t,:,:] == l))
        #filtering overlap in time (TRY)
        """
        if t > 1:
          boolarr1 = arr[t,:,:] == l
          boolarr2 = arr[t-1,:,:] == l
          boolarr = boolarr1*boolarr2
          n = np.count_nonzero(boolarr)
          n_ex = np.count_nonzero(boolarr2)
          if n < n_ex/2:
            ncomponents += 1
            for x in range(t,t+20):
              try:
                boolarr = arr[x,:,:] == l
                arr[x,:,:] = xr.where(boolarr,ncomponents,arr[x,:,:])
              except:
                break
        """
    arr = OrderIndexes(arr)
    #initialize coords for new .nc
    times = pIB_boolean.coords["time"].values
    plev = pIB_boolean.coords["plev"].values
    lon = pIB_boolean.coords["lon"].values
    lat = pIB_boolean.coords["lat"].values
    #initialize dataset object for the new .nc
    pIB_tracked = xr.DataArray(0,coords=[times,plev,lat,lon],dims = pIB_boolean.dims)
    pIB_tracked[:,:,:,:] = 0
    pIB_tracked[:,0,:,:] = arr
    #assign dataset to self.dataset which is now updated
    self.dataset = self.dataset.assign(pIB_tracked = pIB_tracked)

    #output data
    if data_return == False:
      self.dataset.to_netcdf(fn_out)
    if data_return == True:
      return self.dataset
    else:
      return 0
