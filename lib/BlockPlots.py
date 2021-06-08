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
import lib.BlockTools as BlockTools

#plot dimensions in inches
plt.rcParams["figure.figsize"] = (10,10)
#np.set_printoptions(precision=2,threshold=np.inf)

"""
Some data controls for zg and pIB variables.
Set up for following plots.
"""
def check_freq(ds,freq=""):
  try:
    freq = ds[freq]
    string = "data read successfully"
    print(string)
  except:
    string = string = "Error code 2: frequency was not found. " +\
    "Hint: use read_main() to load data with right attributes."
    print(ds)
    print(string)
    return 2
  #preparing horizontal positions
  lats = ds.lat.values
  lons = ds.lon.values
  #adding an element to longitude direction as needed by cartopy stereographic transforms
  if lons[len(lons)-1] != lons[0]:
    lons = np.append(lons,180-10**(-6))#slightly different for avoiding contour bug
    freq = cutil.add_cyclic_point(freq)
  return lons,lats,freq


def check_zg_pIB(ds,zg="",pIB=""):
  try:
    zg = ds[zg]
    pIB = ds[pIB]
    string = "data read successfully"
  except:
    string = "Error code 2: zg variable and pIB were not found. " +\
    "Hint: use read_main() to load data with right attributes."
    print(ds)
    print(string)
    return 2

  #____CHECK GEOP HEIGHT____
  #ERA5 dataset uses geopotential, not height
  if zg.values[0,0,0,0] > 10000:
    	zg = zg/9.80665

  #creating arrays
  dat1 = zg.loc[:,5e+04,:,:].values
  dat2 = pIB.loc[:,5e+04,:,:].values
  #preparing horizontal positions
  lats = ds.lat.values
  lons = ds.lon.values

  #adding an element to longitude direction as needed by cartopy stereographic transforms
  if lons[len(lons)-1] != lons[0]:
    lons = np.append(lons,180-10**(-6))#slightly different for avoiding contour bug
    dat1 = cutil.add_cyclic_point(dat1)
    dat2 = cutil.add_cyclic_point(dat2)
  return lons,lats,dat1,dat2


"""
PlotEventZ500 function produces a plot for a specific day and time of Z500 field
and highlight where the Z500 anomaly index detected a blocking event.
All requested variables have defaults except from output.
"""

def PlotEventZ500(ds,
                  output,
                  day="1979-01-01T09:00:00.000000000",
                  pIB="pIB_boolean", #name of IB matrix variable, may vary from file to file
                  zg="zg", #name of zg500 matrix, may vary from file to file
                  mapcrs = ccrs.AzimuthalEquidistant(central_latitude = 90),
                  datacrs = ccrs.PlateCarree(),
                  extent =[-180,180,20,70],
                  plot_title=""):

  lons,lats,dat1,dat2 = check_zg_pIB(ds,zg,pIB)

  #selecting data from the desired day
  index = BlockTools.GetIndex(ds,"time",day)
  dat1 = dat1[index,:,:]
  dat2 = dat2[index,:,:]

  #multiple plots for having both contours and heatmap
  ax = plt.subplot(111, projection=mapcrs, aspect = 1.0)
  ax.set_extent(extent, datacrs)
  ax.coastlines()

  #define ranges
  cs_range = np.arange(4500,6500,40) #usually it works
  cb_range = [0.,1.]
  print(cb_range)

  #plot contour
  cs = ax.contourf(lons, lats,dat1, cs_range,cmap="jet", transform=datacrs) #first ploot (float)
  cs_b = ax.contour(lons, lats,dat1, cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#second plot for contours
  plt.clabel(cs_b, colors ="black", fontsize = 5, inline ="false")
  cb = ax.contour(lons,lats,dat2,cb_range,cmap = "flag", transform=datacrs)#blocking plot

  # Make some nice titles for the plot
  plt.title(plot_title,fontsize = 12, loc='left')
  plt.colorbar(cs, orientation='horizontal', pad=0, label = "geopotential height" , aspect=100)
  #Export image
  try:
    plt.savefig(output,bbox_inches='tight',dpi=250)
    return 0
  except:
    print("Error code 1: no output file was given or something went wrong with matplotlib")
    return 1

def PlotFreqZ500(ds,\
                 output,\
                 mapcrs=ccrs.PlateCarree(),
                 freq = "TM90_freq",
                 plot_title = "",
                 cmap = plt.cm.BuPu):
  lons,lats,freq_z=check_freq(ds,freq)
  tmp = freq
  freq = freq_z[0,:,:]
  # Set up the projection that will be used for plotting the map
  mapcrs = mapcrs
  # Set up the projection of the data; if lat/lon then PlateCarree is what you want
  datacrs = ccrs.PlateCarree()
  #plot map subplot
  ax = plt.subplot(111, projection=mapcrs)
  ax.set_extent([-180, +180, 0, 90], ccrs.PlateCarree())
  ax.coastlines()
  if tmp == "freq_anomaly":
    cs_range = np.linspace(-7,+7,14,endpoint = True)
    cs = ax.contourf(lons,lats,freq,cs_range,cmap=cmap,transform=datacrs)
  #plot contour
  else:
    cs = ax.contourf(lons, lats, freq,cmap=cmap,transform=datacrs)
  # Make some nice titles for the plot (one right, one left)
  plt.title(plot_title, loc='left')
  plt.colorbar(cs, orientation='horizontal', pad=0, label = "% over the period" , aspect=50)
  #plt.colorbar.label("% over the period")
  #Export image
  try:
    plt.savefig(output,bbox_inches='tight',dpi=250)
    return 0
  except:
    print("Error code 1: no output file was given or something went wrong with matplotlib")
    return 1
"""
Plot Z500 composite over blocking event in a grid point.
If compare_anomaly = True the function expect another variable to compare, for example
temperature or precipitation
"""
def PlotCompositeZ500(ds,
                      add_ds,
                      output,
                      pIB="pIB_boolean",
                      zg="zg",
                      point_coords = [0,0], #coordinates for the point where composite is computed\
                      #latitude -longitude
                      additional = "",
                      mapcrs = ccrs.AzimuthalEquidistant(central_latitude = 90),
                      datacrs = ccrs.PlateCarree(),
                      colorbar_label = "",
                      extent =[-180,180,20,70],
                      plot_title=""):
  #dat1 is zg dat2 is pIB_boolean
  lons,lats,dat1,dat2 = check_zg_pIB(ds,zg,pIB)
  if additional != "":
    try:
      dat3 = add_ds[additional].loc[:,5e+04,:,:].values
      dat3 = cutil.add_cyclic_point(dat3)
      meanT = dat3.mean(0)
    except:
      try:
        dat3 = add_ds[additional].values
        if additional == "ta":
          dat3 = dat3*24*1000
        dat3 = cutil.add_cyclic_point(dat3)
        meanT = dat3.mean(0)
      except:
        print("Error code 1: invalid additional dataset.")
        return 1
    print(np.unique(dat3)[0])

  counter = 0
  ilat= BlockTools.GetIndex(ds,"lat",str(point_coords[0]))
  ilon= BlockTools.GetIndex(ds,"lon",str(point_coords[1]))
  #computing composite
  for t in np.arange(0,len(ds.time.values)):
    if dat2[t,ilat,ilon] == 0:
      dat1[t,:,:] = 0
      if additional != "":
        dat3[t,:,:] = -1.
      counter += 1
  #update dat1
  dat1 = np.ma.masked_equal(dat1,0)#mask array1 when ==0
  dat1 = dat1.mean(0)
  #Compute percentage
  perc = (1-(counter/ds.time.values.shape[0]))*100


  #Initialize group of subplot and draw coastlines
  ax = plt.subplot(111, projection=mapcrs)
  ax.set_extent(extent, ccrs.PlateCarree())
  ax.coastlines()
  if(str(mapcrs) == "ccrs.PlateCarree()"):
     gl = ax.gridlines(crs=datacrs, draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
  #plot contour
  c1_range = np.arange(4500,6500,50)
  c1_b = ax.contour(lons,lats,dat1, c1_range,colors="black",linewidths = 0.1 ,transform=datacrs)#second plot for contours
  plt.clabel(c1_b, colors ="black", fontsize ="smaller", inline ="false")
  #update dat3 if present
  if additional != "":
    dat3 = np.ma.masked_equal(dat3,-1.)
    dat3 = dat3.mean(0) - meanT #temp anomaly
    max = np.amax(np.absolute(dat3))
    c3_range = np.linspace(-max,+max,20,endpoint = True)
    c3 = ax.contourf(lons, lats,dat3,c3_range,cmap="coolwarm", transform=datacrs) #first ploot (float)
    cbar = plt.colorbar(c3, orientation='horizontal', pad=0, label = colorbar_label ,ticks = c3_range[1:-1], aspect=40)
    cbar.ax.tick_params(labelsize=7)
  # Make some nice titles for the plot (one right, one left)
  plt.title(plot_title, loc='left')
  plt.title("% events over the period: " + str(perc)[0:5], loc="right")
  plt.plot(point_coords[1],point_coords[0],"or",transform=datacrs)

  # Export image
  try:
    plt.savefig(output,bbox_inches='tight',dpi=250)
    return 0
  except:
    print("Error code 1: no output file was given or something went wrong with matplotlib")
    return 1

"""
Plot Atmospheric Blocking Tracking.
This function takes a .nc file and a 15 days range. The .nc file must contain an attribute
pIB_tracked which is a time,lon,lat matrix which is zero when  there is no blocking and
n when the nth blocking event is occuring.
The function performs a plot of the 15 days geopotential heights and colors with different
colors the many blocking events.
"""
def PlotTracking(ds,
                 output,
                 starting_day="1997-01-01T09:00:00.000000000",
                 pIB_tracked="pIB_tracked", #name of IB matrix variable, may vary from file to file
                 zg="zg", #name of zg500 matrix, may vary from file to file
                 mapcrs = ccrs.PlateCarree(),
                 datacrs = ccrs.PlateCarree(),
                 extent =[-180,180,0,90]):

  lons,lats,dat1,dat2 = check_zg_pIB(ds,"zg",pIB_tracked)
  #selecting data from the desired days
  index1 = BlockTools.GetIndex(ds,"time",starting_day)
  index2 = index1 + 15
  print(index1,index2)
  dat1 = dat1[index1:index2,:,:]
  dat2 = dat2[index1:index2,:,:]
  #multiple plots for having both contours and flags over different days
  dat2 = BlockTools.OrderIndexes(dat2)
  min = np.amin(dat2) #store min and max for following calculation
  max = np.amax(dat2)
  dat2 = np.ma.masked_equal(dat2,0) #mask array1 when ==0
  #print(dat2)
  # define the colormap
  list = ['w','g','r','m','c','y','k','tab:blue','tab:orange','tab:green',\
          'tab:purple','tab:brown','tab:pink','tab:gray','olive','cyan','gold',\
          'lightcoral','deeppink','lightsteelblue','lime','indigo','palegoldenrod','b']
  N = max+1
  cmap = c.ListedColormap(list[:N])
  # define the bins and normalize
  if N > 1:
    bounds = np.linspace(0,N,N+1)
  else:
    bounds = [0,1,2]
  norm = mpl.colors.BoundaryNorm(bounds, N)
  for i in range(0,15):
    ax = plt.subplot(5,3,i+1, projection=mapcrs, aspect = 1.0)
    ax.set_extent(extent, datacrs)
    ax.coastlines()
    #define ranges
    cs_range = np.arange(4500,6500,80) #usually it works
    #plot contour
    #dat2 is normalized in 0,1
    cs_b = ax.contour(lons, lats,dat1[i,:,:], cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#sec>    plt.clabel(cs_b, colors ="black", fontsize = 5, inline ="false")
    cb = ax.pcolor(lons,lats,dat2[i,:,:],cmap=cmap,vmin = 0,vmax = N,transform=datacrs)#blocking plot
    # create the colorbar
    cb = plt.colorbar(cb, spacing='proportional',orientation='horizontal',ticks=bounds)
    cb.set_label('blocking event label')
    # Make some nice titles for the plot
    plt.title(starting_day[0:8] + starting_day[8:10] + ' + ' + str(i) + ' days',fontsize = 12, loc='left')
  #Export image
  try:
    plt.savefig(output,bbox_inches='tight',dpi=250)
    return 0
  except:
    print("Error code 1: no output file was given or something went wrong with matplotlib")
    return 1
