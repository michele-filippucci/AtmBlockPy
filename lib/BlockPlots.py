"""
______________________________________
//////////////        \\\\\\\\\\\\\\\\
////////                     \\\\\\\\\
||||||||  BLOCKPLOTS LIBRARY  ||||||||
\\\\\\\\_____          ______/////////
\\\\\\\\\\\\\\________////////////////

Author: Michele Filippucci on intership at ISAC-CNR (TO)

This a library complementary to BlockTools.

BlockTools is a set of tools for the analysis of atmospheric blocking in the northern hemisphere.
The index used for atm blocking diagnostic is described in "Davini et al. - 2012 - Bidimensional diagnostics, variability, and 
trends of northern hemisphere blocking". Some differences and features are added: the persistence and area criteria
are applied at the level of tracking. Tracking also allow the user to perform lagrangian analysis.

The function in this script are suitable for quick plotting of typical field used in the analysis of atmospheric blocking. 

The requirements for this library are:
Python              3.8.10
xarray              0.18.2
numpy               1.20.3
scipy               1.6.3
tqdm                4.61.1
matplotlib          3.4.2
cartopy             0.19.0.post1

______________________________________
//////////////        \\\\\\\\\\\\\\\\
////////                     \\\\\\\\\
||||||||  LIST OF FUNCTIONS:  ||||||||
\\\\\\\\_____          ______/////////
\\\\\\\\\\\\\\________////////////////

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ bp_automatic_colorbar(symmetric = False,starting_val = 0,ending_val = 0,n_bins = 10) _ _ _ 

output: colors_bins (list),ticks (list)

This function produces an automatic colorbar that can be used with every automatic color scheme offered by matplotlib (es: "Blues", "Reds",
"jet", "bwr").
symmetric: useful for anomaly plot where "bwr" is used and we want white to refer to zero values
starting_val: first element of colorbar and ticks
ending_val: ...
n_bins: ...

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ bp_colorbar(phys_quantity,half = False,symmetric = False,starting_val = -10,ending_val = +10, n_bins = 10) _ _ _ 

This function is actually a list of suitable colorbar for a range of possible plots. phys_quantities supported can be found in the function
implementation.
The other arguments are the same as bp_automatic_colorbar in case the phys_quantity is not implemented.

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ PlotField(ds,output = "",field="zg",mapcrs = ccrs.AzimuthalEquidistant(central_latitude = 90),
                  datacrs = ccrs.PlateCarree(),
                  extent =[-180,180,20,70],
                  box_vertices = [[0,0],[0,0],[0,0],[0,0]],
                  subplot_grid = 11,
                  subplot_idx = 1,
                  plot_title=""): _ _ _ 
               
This function produces a plot of a field on a map.
ds : input dataset
output: plot destination
field: name of the field contained in the dataset
mapcrs: projection of the map
datacrs: projection of the data
extent: extension in degrees of the plot
box_vertices: in case it's a composite it is possible to plot the box where the composite was computed.
subplot_grid : if different from 11 the function returns an axis object instead of a saved plot. This way it is easy to produce
multiplots
plot_title: ...

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

_ _ _ PlotComposite(  ds,
                      output,
                      add_ds= "",
                      zg_mean="zg_mean",
                      var_anomaly = "var_anomaly",
                      cb_endpoint = 0,
                      point_coords = [0,0], #coordinates for the point where composite is computed\
                      #latitude -longitude
                      additional = "", #additional quantity : for example temperature. must be stored in add_ds
                      mapcrs = ccrs.AzimuthalEquidistant(central_latitude = 90),
                      datacrs = ccrs.PlateCarree(),
                      colorbar_label = "",
                      extent =[-180,180,20,70],
                      plot_title="") _ _ _
                      
This function is able to plot a composite. This is a complex task and this function is hardly useful.
The arguments are pretty similar to those of PlotField

_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

PlotTracking(    ds,
                 output,
                 starting_day="1997-01-01T09:00:00.000000000",
                 pIB_tracked="pIB_tracked", #name of IB matrix variable, may vary from file to file
                 zg="zg", #name of zg500 matrix, may vary from file to file
                 mapcrs = ccrs.PlateCarree(),
                 datacrs = ccrs.PlateCarree(),
                 extent =[-180,180,0,90]):
                 
This function is able to plot a series of 15 days and highliting the various events taking place during each day.
                 
_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/-\_/

PlotLabel(ds,
              output,
              label,
              shuffle = False,
              pIB_tracked="pIB_tracked", #name of IB matrix variable, may vary from file to file
              zg="zg", #name of zg500 matrix, may vary from file to file
              mapcrs = ccrs.PlateCarree(),
              datacrs = ccrs.PlateCarree(),
              extent =[-180,180,0,90])
              
This function is able to produce a plot of the first 15 days of life of an event with an arbitrary label. Arguments are self-explanatory.
shuffle : if shuffle = True the plot is shifted of 15 days.


"""

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

#my own classes
import sys
sys.path.append("/home/guest/work/michele/prog/GitPROG/AtmBlockPy/lib")
import BlockTools as BlockTools


#plot dimensions in inches
plt.rcParams["figure.figsize"] = (9,15)
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
  try:
    dat1 = zg.loc[:,5e+04,:,:].values
    dat2 = pIB.loc[:,5e+04,:,:].values
  except:
    dat1 = zg.values[:,0,:,:]
    dat2 = pIB.values[:,0,:,:]
  #preparing horizontal positions
  lats = ds.lat.values
  lons = ds.lon.values

  #adding an element to longitude direction as needed by cartopy stereographic transforms
  if lons[len(lons)-1] != lons[0]:
    lons = np.append(lons,180-10**(-6))#slightly different for avoiding contour bug
    dat1 = cutil.add_cyclic_point(dat1)
    dat2 = cutil.add_cyclic_point(dat2)
  return lons,lats,dat1,dat2

def check_field(ds,field=""):
  try:
    field_value = ds[field]
    string = "data read successfully"
  except:
    string = "Error code 2: field variable was not found. " +\
    "Hint: use read_main() to load data with right attributes."
    print(ds)
    print(string)
    return 2

  #____CHECK GEOP HEIGHT____
  #ERA5 dataset uses geopotential, not height
  if field == "zg" and (field_value.values > 10000).any():
    field_value = field_value/9.80665

  #preparing horizontal positions
  lats = ds.lat.values
  lons = ds.lon.values
  dat1 = field_value.values

  #adding an element to longitude direction as needed by cartopy stereographic transforms
  if lons[len(lons)-1] != lons[0]:
    lons = np.append(lons,180-10**(-6))#slightly different for avoiding contour bug
    dat1 = cutil.add_cyclic_point(dat1)
  return lons,lats,dat1

"""
A function aimed to produce a suitable colorbar for a plot.
It returns two ranges: the first one containing the class limits and the second containing the tick
parameter.
"""
def bp_automatic_colorbar(symmetric = False,
                starting_val = 0, #useless if colorbar is symmetric
                ending_val = 0,  
                n_bins = 10
                ):
  if symmetric:
    starting_val = - ending_val
    
  if symmetric and n_bins%2 != 0:
    n_bins += 1
    print("n_bins must be even so that colors are odd")
 
  colors_bins = np.linspace(starting_val,ending_val,n_bins,endpoint = True)
  
  #the following complicated loop provides a symmetric colorbar. Every q color bin a tick is plotted and when zero is crossed
  #the counting starts from the right side.
  
  tmp=0
 
  if symmetric:
    
    ticks_idx = []
    count = 0
    if n_bins < 12:
      q = 2
    if n_bins >= 12 and n_bins < 24:
      q = 3
    if n_bins >= 24:
      q = 4
    for val in colors_bins:
      if tmp*val < 0:
        count = 0
      if count%q == 0:
        if val < 0:
          ticks_idx.append(count)
        if val > 0:
          ticks_idx.append(len(colors_bins)-count-1)
      
      tmp = val.copy()
      count +=1
    ticks = colors_bins[ticks_idx]
    
  #we want 8 ticks
  #n_colors_bins/x = n_ticks --> x = int(n_colors_bins/n_ticks)
  if not symmetric:
    n_ticks = 8
    ticks = colors_bins[range(0,len(colors_bins),int(len(colors_bins)/n_ticks))]
    
  #DEBUG
  #print(colors_bins,ticks)
  
  return colors_bins,ticks

def  bp_colorbar(phys_quantity,
                 half = False,
                symmetric = False,
                starting_val = -10, #useless if colorbar is symmetric
                ending_val = +10,  
                n_bins = 10
                ):
  cb = []
  ticks = []
  cb_label = ""
  cmap = ""
  if phys_quantity not in ("m_heat_flux","kin_en","eady_growth_rate"):
    cb,ticks = bp_automatic_colorbar(symmetric,starting_val,ending_val,n_bins)
    cb_label = "unknown"
    cmap = "Reds"
    
  """
  TA
  """
  
       
  if phys_quantity == "ta":
    cb = np.arange(250,350,5) #500
    if half:
      cb = np.arange(250,350,10) #500
    ticks = np.arange(250,350,10) #500
    cb_label = "T (K)"
    cmap = "Reds"
  if phys_quantity == "ta_an":
    cb = np.arange(-7,+8,2)
    if half:
      cb = np.arange(-7.5,+8,2) #500
    ticks = np.arange(-5,+6,2) #500
    cb_label = "T (K)"
    cmap = "bwr"
   
  
  """
  M_HEAT_FLUX
  """
  
  if phys_quantity == "m_heat_flux":
    cb = np.arange(-1,17,1)
    if half:
      cb = np.arange(-1,17,2)
    ticks = np.arange(-1,17,2)
    cb_label = "(m/s)*K"
    cmap = "Reds"
  if phys_quantity == "m_heat_flux_an":
    cb = np.arange(-1.3,1.5,0.2)
    if half:
      cb = np.arange(-1.3,1.5,0.4)
    ticks = np.arange(-1.2,1.4,0.4)
    cb_label = "(m/s)*K"
    cmap = "bwr"
    
  """
  KIN_EN
  """
    
  if phys_quantity == "kin_en":
    cb = np.arange(0,330,30) #250
    if half:
      cb = np.arange(0,330,60) #250
    ticks = np.arange(0,330,60) #250
    cb_label = "(m/s)^2"
    cmap = "Reds"
  if phys_quantity == "kin_en_an":
    #250hPa
    cb = np.arange(-19.5,20,3)
    if half:
      cb = np.arange(-19.5,20,6) #250
    ticks = np.arange(-19.5,20,3) #250
    cb_label = "(m/s)^2"
    cmap = "bwr"
    
  """
  UA
  """
    
  if phys_quantity == "ua":
    #cb = np.arange(-14,15,4) #850
    cb = np.arange(-70,71,20) #250
    if half:
      #cb = np.arange(-14,15,4) #850
      cb = np.arange(-70,71,20) #250
    #ticks = np.arange(-14,15,4) #850
    ticks = np.arange(-70,71,20) #250
    cb_label = "m/s"
    cmap = "bwr"
  if phys_quantity == "ua_an":
    #850hPa
    #cb = np.arange(-2.1,2.2,0.20)
    #250hPa
    cb = np.arange(-4.2,4.3,0.40)
    if half:
      #cb = np.arange(-2.1,2.2,0.60)#850
      cb = np.arange(-4.2,4.3,1.2) #250
    #ticks = np.arange(-2.1,2.2,0.60) #850
    ticks = np.arange(-4.2,4.3,1.2) #250
    cb_label = "m/s"
    cmap = "bwr"
    
  """
  EADY GROWTH RATE
  """
  
       
  if phys_quantity == "eady_growth_rate":
    cb = np.arange(0,1.2,0.1) #850
    if half:
      cb = np.arange(0,1.2,0.2) #850
    ticks = np.arange(0,1.2,0.2	) #850
    cb_label = "(1/day)"
    cmap = "Reds"
  if phys_quantity == "eady_growth_rate_an":
    #250hPa
    cb = np.arange(-0.1,0.1,0.02)
    if half:
      cb = np.arange(-0.1,0.1,0.02) #850
    ticks = np.arange(-0.1,0.1,0.02) #850
    cb_label = "(1/day)"
    cmap = "bwr"
   
  
  """
  FREQUENCY
  """
  
       
  if phys_quantity == "freq":
    cb = np.arange(0,25,2) #850
    if half:
      cb = np.arange(0,25,4) #850
    ticks = np.arange(0,25,4) #850
    cb_label = "(% days)"
    cmap = "Purples"
  if phys_quantity == "freq_an":
    #250hPa
    cb = np.arange(-9.5,9.6,1)
    #cb = np.arange(-3.25,3.3,0.5)
    if half:
      cb = np.arange(-9.75,9.8,2) #850
      #cb = np.arange(-3.25,3.3,1)
    ticks = (-9.5,-7.5,-5.5,-3.5,-1.5,1.5,3.5,5.5,7.5,9.5)
    #ticks = (-3.25,-2.25,-1.25,0,+1.25,+2.25,+3.25)
    cb_label = "(% days)"
    cmap = "bwr"
    
  """
  ZG
  """
  
       
  if phys_quantity == "zg":
    cb = np.arange(5500,6200,40) #500
    if half:
      cb = np.arange(5500,6200,80) #500
    ticks = np.arange(5500,6200,160) #500
    cb_label = "geopotential height (m)"
    cmap = "jet"
  if phys_quantity == "zg_an":
    cb = np.arange(-36.25,+36.5,2.5)
    if half:
      cb = np.arange(-57.5,+58,10) #500
    ticks = (-32.5,-22.5,-12.5,-2.5,2.5,12.5,22.5,32.5) #500
    cb_label = "geopotential height (m)"
    cmap = "bwr"
   
  
  
  return cb,ticks,cb_label,cmap
  
"""
PlotEventZ500 function produces a plot for a specific day and time of Z500 field
and highlight where the Z500 anomaly index detected a blocking event.
All requested variables have defaults except from output.
"""

def PlotField(    ds,
                  output = "",
                  field="zg", #name of field500 matrix, may vary from file to file
                  mapcrs = ccrs.AzimuthalEquidistant(central_latitude = 90),
                  datacrs = ccrs.PlateCarree(),
                  extent =[-180,180,20,70],
                  box_vertices = [[0,0],[0,0],[0,0],[0,0]],
                  subplot_grid = 11,
                  subplot_idx = 1,
                  plot_title=""):

  lons,lats,dat1 = check_field(ds,field)

  #multiple plots for having both contours and heatmap
  ax = plt.subplot(int(str(subplot_grid)+str(subplot_idx)), projection=mapcrs, aspect = 1.0)
  ax.set_extent(extent, datacrs)
  ax.coastlines()
  gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')

  #define colorbar
  cs_range,ticks,cb_label,cmap = bp_colorbar(field)
  
  #plot contour
  cs = ax.contourf(lons, lats,dat1, cs_range,cmap=cmap, transform=datacrs) #first ploot (float)
  cs_b = ax.contour(lons, lats,dat1, cs_range,colors="gray",linewidths = 0.1 ,transform=datacrs)#second plot for contours
  #plt.clabel(cs_b, colors ="black", fontsize = 5, inline ="false")
  
  #plot box if needed
  if box_vertices[0] != box_vertices[1]:
    polygon = plt.Polygon(box_vertices,fill = False,color="Green",transform=datacrs)
    ax.add_patch(polygon)
    ax.patch.set_alpha(1)
    

  # Make some nice titles for the plot
  plt.title(plot_title,fontsize = 12, loc='left')
  plt.colorbar(cs, orientation='horizontal',ticks=ticks,pad=0.08, label = cb_label , aspect=40)
  #Export image
  if subplot_grid != 11:
    return ax
  try:
    plt.savefig(output,bbox_inches='tight',dpi=250)
    return 0
  except:
    print("Error code 1: no output file was given or something went wrong with matplotlib")
    return 1

"""
Plot composite over blocking event in a grid point.
The function need a dataset containing zg and a second dataset containing the second variable
"""
def PlotComposite(ds,
                      output,
                      add_ds= "",
                      zg_mean="zg_mean",
                      var_anomaly = "var_anomaly",
                      cb_endpoint = 0,
                      point_coords = [0,0], #coordinates for the point where composite is computed\
                      #latitude -longitude
                      additional = "", #additional quantity : for example temperature
                      mapcrs = ccrs.AzimuthalEquidistant(central_latitude = 90),
                      datacrs = ccrs.PlateCarree(),
                      extent =[-180,180,20,70],
                      plot_title=""):
  #dat1 is zg dat2 is pIB_boolean
  lons,lats,dat1 = check_field(ds,zg_mean)
  if additional != "":
    try:
      dat3 = add_ds[additional].loc[:,5e+04,:,:].values
      dat3 = cutil.add_cyclic_point(dat3)
    except:
      try:
        dat3 = add_ds[additional].values
        if additional == "ta":
          dat3 = dat3*24*1000
        dat3 = cutil.add_cyclic_point(dat3)
      except:
        print("Error code 1: invalid additional dataset.")
        return 1
    #print(np.unique(dat3)[0])

  #Initialize group of subplot and draw coastlines
  ax = plt.subplot(111, projection=mapcrs)
  ax.set_extent(extent, ccrs.PlateCarree())
  ax.coastlines()
  gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
  if(str(mapcrs) == "ccrs.PlateCarree()"):
     gl = ax.gridlines(crs=datacrs, draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
  #plot contour
  c1_range,ticks,cb_label,cmap = bp_colorbar(zg_mean) #usually it works
  c1_b = ax.contour(lons,lats,dat1, c1_range,colors="black",linewidths = 0.1 ,transform=datacrs)#second plot for contours
  plt.clabel(c1_b, colors ="black", fontsize ="smaller", inline ="false")
  
  #plot dat3 if present
  if additional != "":
    if cb_endpoint != 0:
      max = cb_endpoint
    if cb_endpoint == 0:
      max = np.amax(np.absolute(dat3))
    c3_range,ticks,cb_label,cmap=bp_colorbar(additional) #usually it works
    c3 = ax.contourf(lons, lats,dat3,c3_range,cmap="bwr", transform=datacrs) #first ploot (float)
    cbar = plt.colorbar(c3, orientation='horizontal',pad=0.08, label = cb_label ,ticks = ticks, aspect=60)
    cbar.ax.tick_params(labelsize=7)
    
  # Make some nice titles for the plot (one right, one left)
  plt.title(plot_title, loc='left')
  plt.plot(point_coords[1],point_coords[0],"or",transform=datacrs)

  # Export image
  try:
    plt.savefig(output,bbox_inches='tight')
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
    #cb = plt.colorbar(cb, spacing='proportional',orientation='horizontal',ticks=bounds)
    #cb.set_label('blocking event label')
    # Make some nice titles for the plot
    plt.title(starting_day[0:8] + starting_day[8:10] + ' + ' + str(i) + ' days',fontsize = 12, loc='left')
  #Export image
  try:
    plt.savefig(output,bbox_inches='tight',dpi=250)
    return 0
  except:
    print("Error code 1: no output file was given or something went wrong with matplotlib")
    return 1
  
def PlotLabel(ds,
              output,
              label,
              shuffle = False,
              pIB_tracked="pIB_tracked", #name of IB matrix variable, may vary from file to file
              zg="zg", #name of zg500 matrix, may vary from file to file
              mapcrs = ccrs.PlateCarree(),
              datacrs = ccrs.PlateCarree(),
              extent =[-180,180,0,90]):
  lons,lats,dat1,dat2 = check_zg_pIB(ds,"zg",pIB_tracked)
  #selecting data from the desired days
  #create a list of indexes
  for t in range(len(ds["time"].values)):
    if label in dat2[t,:,:]:
      index1 = t
      break
  time = ds["time"].values
  starting_day = str(time[index1])
  if shuffle == True:
    index1 = index1 + 15
  index2 = index1 + 15
  dat1 = dat1[index1:index2,:,:]
  dat2 = dat2[index1:index2,:,:]
  print(index1,index2)
  bl = dat2!=label
  dat2[bl] = 0
  #multiple plots for having both contours and flags over different days
  min = np.amin(dat2) #store min and max for following calculation
  max = np.amax(dat2)
  dat2 = np.ma.masked_equal(dat2,0) #mask array1 when ==0
  for i in range(0,15):
    ax = plt.subplot(5,3,i+1, projection=mapcrs, aspect = 1.0)
    ax.set_extent(extent, datacrs)
    ax.coastlines()
    #define ranges
    cs_range = np.arange(4500,6500,80) #usually it works
    N=1
    #plot contour
    #dat2 is normalized in 0,1
    cs_b = ax.contour(lons, lats,dat1[i,:,:], cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#sec>    plt.clabel(cs_b, colors ="black", fontsize = 5, inline ="false")
    cb = ax.pcolor(lons,lats,dat2[i,:,:],cmap="Oranges",vmin = 0,vmax = N,transform=datacrs)#blocking plot
    # create the colorbar
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

  
