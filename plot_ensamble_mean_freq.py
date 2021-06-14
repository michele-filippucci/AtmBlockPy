import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import colors as c
import matplotlib as mpl
import numpy as np
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#my own class
from lib.BlockTools import BlockTools
import lib.BlockTools as BT
import lib.BlockPlots as BP

compute_freq = True
MGF = False
version = "las"

if MGF == True:
  mgf = "_mgf_"
else:
  mgf = ""

output = "/home/guest/work/michele/prog/plots/ensamble/ens"
reanalysis = "/home/guest/work/michele/data/ERA5/Z500/"+\
     "ERA5_Z500_day_djfm_r144x73_500hPa_northem_1979-2008.nc"
fn_out = "/home/guest/work/michele/data/ERA5/processed/"+\
         "ERA5_pIB_daily_djfm_northem_1979-2008.nc"
climatology="/home/guest/work/michele/data/ERA5/processed/"+\
         mgf+version + "_TM90_clim_djfm_northem_1979-2008.nc"
#THIS PART PRODUCES THE DATA
if compute_freq == True:
  var = []
  for ens_number in range(10):
    print("Starting iteration N: " + str(ens_number))
    string = "/home/guest/work/michele/data/SPHINX/"+version+str(ens_number)+"/zg/total_northem_djfm_"+version+str(ens_number)+".nc"
    print("Opening file:" + string)
    ds = xr.load_dataset(string)
    bloc = BlockTools()
    bloc.load_data(ds)
    ds = bloc.TM90(data_return=True,freq_also=True,mer_gradient_filter=MGF)
  #  BP.PlotFreqZ500(ds=ds,output=output+str(ens_number),plot_title=output+str(ens_number))
    freq = ds["TM90_freq"].values
    if ens_number == 0:
      tuple = np.zeros(np.shape(freq))
    var.append(freq)
    tuple += freq

  #initialize coordinates for dataset
  zg = ds["zg"]
  plev = zg.coords["plev"].values
  lon = zg.coords["lon"].values
  lat = zg.coords["lat"].values

  #produce ERA5 data
  pIB = BlockTools()
  pIB.read(reanalysis)
  ds2 = pIB.TM90(data_return=True,freq_also = True,mer_gradient_filter = MGF)

  #produce dataset
  freq_ensamble = xr.DataArray(0,coords=[plev,lat,lon],dims = ["plev","lat","lon"])
  freq_ensamble[:,:,:] = tuple/10

  freq_anomaly = xr.DataArray(0,coords=[plev,lat,lon],dims = ["plev","lat","lon"])
  freq_anomaly[:,:,:] = freq_ensamble.values-ds2["TM90_freq"].values

  freq_ERA5 = xr.DataArray(0,coords=[plev,lat,lon],dims = ["plev","lat","lon"])
  freq_ERA5[:,:,:] = ds2["TM90_freq"].values

  sm = np.zeros(np.shape(var[0]))
  for i in range(10):
    var[i] = (var[i]-freq_ensamble)**2
    sm += var[i]
  variance = sm/10

  freq_variance = xr.DataArray(variance,coords=[plev,lat,lon],dims = ["plev","lat","lon"])
  freq_variance[:,:,:] = variance
  print(np.shape(freq_variance))
  print(np.shape(variance))

  ds = xr.Dataset()

  ds = ds.assign(freq_ensamble=freq_ensamble)
  ds = ds.assign(freq_anomaly=freq_anomaly)
  ds = ds.assign(freq_ERA5=freq_ERA5)
  ds = ds.assign(freq_variance=freq_variance)
  ds.to_netcdf(climatology) #save dataset with climatology
if compute_freq == False:
  ds = xr.load_dataset(climatology)
#THIS PART PRODUCES THE PLOT
#now plot anomaly
BP.PlotFreqZ500(ds=ds,output=output+version+mgf+"ERA5",mapcrs = ccrs.NorthPolarStereo(),freq="freq_ERA5",plot_title="1979-2008 ERA5 : TM90 Freq")
BP.PlotFreqZ500(ds=ds,output=output+version+mgf+"Ensamble",mapcrs = ccrs.NorthPolarStereo(),freq="freq_ensamble",plot_title = "1979/2008 "+version+" : Ensamble Mean TM90 Freq")
BP.PlotFreqZ500(ds=ds,output=output+version+mgf+"Era5Anomaly",mapcrs = ccrs.NorthPolarStereo(),freq="freq_anomaly",plot_title = "1979/2008 "+version+" - ERA5 : TM90 Freq",cmap = plt.cm.bwr)
print(version+" :Integrated absolute anomaly over every grid point:")
print(sum(sum(sum(np.absolute(ds["freq_anomaly"].values)))))
print("Grid points structure:")
print(np.shape(ds["freq_anomaly"].values))
