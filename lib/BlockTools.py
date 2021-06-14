import numpy as np
import xarray as xr
from scipy.ndimage.measurements import label
from scipy.ndimage.measurements import center_of_mass
import cartopy.util as cutil
import sys

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
    arr = xr.where(arr == newarr[i], newval[i],arr)
  return arr

"""
This function returns the coordinates in space of the center of mass
of a blocking event. The object returned is a list of the coordinates
at each time step.
"""

def CenterofMass(tuple,label,grid=2.5):
#  print("CIAOOOOOOOOOOOOOOOOOOO")
  time = np.shape(tuple)[0]
  x = []
  y = []
  bool = xr.where(tuple == label,1.,0.)
  for t in range(time):
    try:
      cm = center_of_mass(bool[t,:,:])
    except:
      print("The label " + label + " isn't present")
      break
#    if len(np.unique(bool[t,:,:]))==1:
#      break
    cm = np.array([cm[0]*2.5,cm[1]*2.5-180])
    x.append(cm[0])
    y.append(cm[1])
  return x,y

class BlockTools(object):
  num_of_BlockTools = 0

  def __init__(self):
    BlockTools.num_of_BlockTools += 1

  def load_data(self,ds):
    self.dataset = ds

  def read(self,filename):
    self.dataset = xr.load_dataset(filename)
    return self.dataset

  """
  This function is capable of creating an nc file identical (located in fn_out)
  to the BlockTools.dataset with additional attributes:
  pIB_boolean and (when freq_also == True) pIB_frequencies
  As an alternative you can change the flag "data_return" and the function
  will return a dataset object from the class xarray containing the same
  additional attributes
  """

  def TM90(self,fn_out = "",\
           data_return = False,\
           freq_also = False,\
           mer_gradient_filter = False,\
           long_filter = False):

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

    #define XArray using the costructor for an appropriate output
    times = zg.coords["time"].values
    plev = zg.coords["plev"].values
    lon = zg.coords["lon"].values
    lat = zg.coords["lat"].values

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
    if mer_gradient_filter == False:
      TuplepIB = xr.where(GHGN < -10.0, 1.0, 0.0) * xr.where(GHGS > 0., 1.0 , 0.0)
    #filter for meridional gradient
    if mer_gradient_filter == True:
      GHGS2 = (+ zg.loc[:,:,15:60,:].values - zg.loc[:,:,0:45,:].values)/15.0
      TuplepIB = xr.where(GHGS2 < -5,1.0,0.0)*xr.where(GHGN < -10.0, 1.0, 0.0)*\
                 xr.where(GHGS > 0., 1.0 , 0.0)
    check = TuplepIB
    #15 degrees continuous longitude filter
    if long_filter == True:
      temp1 = TuplepIB
      temp2 = temp1
      shift = 3
      for i in range(len(lon)):
        #shift 3 to the right
        if i < shift:
          temp1[:,:,:,i] = TuplepIB[:,:,:,len(lon)-shift+i]
        else:
          temp1[:,:,:,i] = TuplepIB[:,:,:,i-shift]
        #shift 3 to the left
        if i < len(lon)-shift:
          temp2[:,:,:,i] = TuplepIB[:,:,:,i+shift]
        else:
          temp2[:,:,:,i] = TuplepIB[:,:,:,i-len(lon)+shift]
      print(temp1[50,0,:,:]*temp2[50,0,:,:]-TuplepIB[50,0,:,:])
      TuplepIB = temp1*temp2
      print(np.unique(TuplepIB-check))
      #pIB = pIB.loc[:,:,:,(-180-7.5):(180-7.5)]*pIB.loc[:,:,:,(-180+7.5):(180+7.5)]

    #define XArray using the costructor for an appropriate output
    pIB = xr.DataArray(0,coords=[times,plev,lat,lon],dims = zg.dims)
    pIB.loc[:,:,30:75,:] = TuplepIB

    self.dataset = self.dataset.assign(TM90 = pIB)
    if freq_also == True:
      pIB_f = sum(pIB)*100/pIB.values.shape[0]
      self.dataset = self.dataset.assign(TM90_freq = pIB_f)
    if data_return == False:
      self.dataset.to_netcdf(fn_out)
    if data_return == True:
      return self.dataset
    else:
      return 0

  """
  Tibaldi and Molteni Index
  This function takes a .nc file containing z500 variable and computes the Tibaldi
  and Monteni index for the latitude 60° N.
  It outputs a .dat file containing the design matrix (features, boolean label) needed
  for training a neural network.
  """

  def TM(self,
         output):
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
        zg = zg.values/9.80665
    #.values gives tuples
    print(self.dataset["lat"])
    N = GetIndex(self.dataset,"lat","75.0") #north
    C = GetIndex(self.dataset,"lat","60.0") #center
    S = GetIndex(self.dataset,"lat","45.0") #south
    print(N,C,S)
    file = open(output, "a")
    for i in range(len(self.dataset["time"])):
      for j in range(len(self.dataset["lon"])):
        string = ""
        for k in range(12):
          string += str(zg[i,0,S+k,j]) + " "
        GHGS = (+ zg[i,0,C,k] - zg[i,0,S,k])/15.0
        GHGN = (- zg[i,0,C,k] + zg[i,0,N,k])/15.0
        flag = int(GHGN < -10.0) * int(GHGS > 0.)
        string += str(flag)
#        print(string)
        file.write(string + "\n")
    file.close()
    return 0


  """
  Contour Tracking
  This function takes a .nc file containing the pIB_boolean attribute from the function
  boolean_pIB and creates a new .nc file containing an additional attribute called 
  pIB_tracked which is zero when blocking is not occuring and is n when the nth blocking
  event is occuring.
  3D version uses label method on a 3D matrix (time,longitude,latitude)
  2D version uses label method in a foor loop (on t) on a 2D matrix (longitude,latitude)
  """
  def ContourTracking3D(self,fn_out = "",var_name = "pIB_boolean",data_return = False):
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
    #applying some filters
    for t in np.arange(0,len(self.dataset.time.values-1)):
      bool = arr[t,:,:] > 0
      list = np.unique(arr[t,bool])
      for l in list:
        boolarr = arr[t,:,:] == l
        n = np.count_nonzero(boolarr)
        #filtering cluster dimension
        if n < 9:
          arr[t,:,:] = xr.where(boolarr, 0,arr[t,:,:])

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

  def ContourTracking2D(self,fn_out = "",var_name = "pIB_boolean",data_return = False, pers = 0):
    if fn_out=="" and data_return==False:
      string = "Specify the kind of output you want"
      print(string)
      return 0
    try:
      pIB_boolean = self.dataset[var_name]
    except:
      print("Error Code 1: dataset not valid. The variable " + var_name + " cannot be found")
      return 1

    #loop over time
    max = 0
    lastlat = len(self.dataset.lat.values) -1
    lastlon = len(self.dataset.lon.values) -1
    for t in np.arange(0,len(self.dataset.time.values-1)):
      if t > 0:
        tmp = np.amax(arr[t-1,:,:])
        if max < tmp: #update maximum value in matrix
          max = tmp
      arr = pIB_boolean.values[:,0,:,:]

      #label method from scipy.ndimage.measurements is used
      structure = [[0,1,0],\
                   [1,1,1],\
                   [0,1,0]] #this matrix defines what is defined as neighbour

      #neighbour points are labeled with the same sequential number
      arr[t,:,:],ncomponents=label(arr[t,:,:],structure)
      arr[t,:,:] = xr.where(arr[t,:,:] > 0, arr[t,:,:] + max , arr[t,:,:])

      #making it periodic in longitude
      for j in range(0,lastlat):
        if arr[t,j,lastlon] > 0 and arr[t,j,0] > 0:
          arr[t,:,:] = xr.where(arr[t,:,:] == arr[t,j,lastlon], arr[t,j,0], arr[t,:,:])

      #applying some filters
      bool = arr[t,:,:] > 0
      comp = np.unique(arr[t,bool])
      for l in comp:
        boolarr = arr[t,:,:] == l
        n = np.count_nonzero(boolarr)

        #filtering cluster dimension
        if n < 9:
          arr[t,:,:] = xr.where(boolarr, 0,arr[t,:,:])

      """
      TRACKING IN TIME
      """
      if t > 0:
        lbl = 1
        bool1 = arr[t-1,:,:] > 0
        bool2 = arr[t,:,:] > 0
        comp1 = np.unique(arr[t-1,bool1])
        comp2 = np.unique(arr[t,bool2])
        for l1 in comp1:
          boolarr1 = arr[t-1,:,:] == l1
          for l2 in comp2:
            boolarr2 = arr[t,:,:] == l2
            boolarr = boolarr1*boolarr2
            n = np.count_nonzero(boolarr)
            n_ex = np.count_nonzero(boolarr1)
            n_new = np.count_nonzero(boolarr2)
            if n > n_ex/2 or n > n_new/2: #50% overlap
              #new label which is always different
              arr[t,:,:] = xr.where(boolarr2,l1,arr[t,:,:])
    arr[:,:,:] = OrderIndexes(arr[:,:,:])

    """
    PERSISTENCE MODULE
    """
    if pers > 0:
      counter = 0
      safe = []
      for t in np.arange(0,len(self.dataset.time.values)):
        bool1 = arr[t,:,:] > 0
        try:
          bool2 = arr[t+pers,:,:] > 0
        except:
          arr[t:,:,:] = 0 #not possible to check so out of output
          print("exited with " + str(len(self.dataset.time.values)-t-pers) + " elements remaining")
          print(str(counter) + " blocking events where found")
          break
        comp1 = np.unique(arr[t,bool1]) #labels at day t
        comp2 = np.unique(arr[t+pers,bool2]) #labels at day t+pers
        for l1 in comp1:
          if not l1 in safe:
            if not l1 in comp2: #if pers days after there is no l1 the event is deleted
              arr[:,:,:] = xr.where(arr[:,:,:]==l1,0,arr[:,:,:])
            else:
              safe.append(l1) #if pers days after there is l1 the event is saved
              counter += 1

    arr = OrderIndexes(arr)
    bool = arr > 0
    print("number of labels: " + str(len(np.unique(arr[bool]))))

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
