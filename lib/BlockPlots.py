import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#my own class
from lib.BlockTools import BlockTools
import lib.BlockTools as BlockTools

#plot dimensions in inches
plt.rcParams["figure.figsize"] = (16,9)
#np.set_printoptions(precision=2,threshold=np.inf)

class BlockPlots(object):
  num_of_BlockPlots = 0
  plot_title = ""
  def __init__(self,plot_title):
    self.plot_title = plot_title
    self.num_of_BlockPlots += 1

  def read_main(self,filename):
    self.main_dataset = xr.load_dataset(filename)

  def read_additional_ds(self,filename):
    self.additional_dataset = xr.load_dataset(filename)

  """
  Some data controls for zg and pIB variables.
  Set up for following plots.
  """
  def check_zg_pIB(self,zg="",pIB=""):
    try:
      zg = self.main_dataset[zg]
      pIB = self.main_dataset[pIB]
      string = "data read successfully"
    except:
      string = "Error code 2: zg variable and pIB were not found. " +\
      "Hint: use read_main() to load data with right attributes."
      print(self.main_dataset)
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
    lats = self.main_dataset.lat.values
    lons = self.main_dataset.lon.values

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

  def PlotEventZ500(self,
           output,
           day="1979-01-01T09:00:00.000000000",
           pIB="pIB_boolean", #name of IB matrix variable, may vary from file to file
           zg="zg", #name of zg500 matrix, may vary from file to file
           mapcrs = ccrs.AzimuthalEquidistant(central_latitude = 90),
           datacrs = ccrs.PlateCarree(),
           extent =[-180,180,20,70]):

    lons,lats,dat1,dat2 = self.check_zg_pIB(zg,pIB)

    #selecting data from the desired day
    index = BlockTools.GetIndex(self.main_dataset,"time",day)
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
    plt.title(self.plot_title,fontsize = 12, loc='left')
    plt.colorbar(cs, orientation='horizontal', pad=0, label = "geopotential height" , aspect=100)
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
  def PlotCompositeZ500(self,
           output,
           pIB="pIB_boolean",
           zg="zg",
           point_coords = [0,0], #coordinates for the point where composite is computed\
                                 #latitude -longitude
           additional = "",
           mapcrs = ccrs.AzimuthalEquidistant(central_latitude = 90),
           datacrs = ccrs.PlateCarree(),
           colorbar_label = "",
           extent =[-180,180,20,70]):

    lons,lats,dat1,dat2 = self.check_zg_pIB(zg,pIB)

    if additional != "":
      try:
        dat3 = self.additional_dataset[additional].loc[:,5e+04,:,:].values
        #print("ok")
        dat3 = cutil.add_cyclic_point(dat3)
        meanT = dat3.mean(0)
      except:
        try:
          dat3 = self.additional_dataset[additional].values*1000
          dat3 = cutil.add_cyclic_point(dat3)
          meanT = dat3.mean(0)
        except:
          print("Error code 1: invalid additional dataset.")
          return 1
      #print(dat3)

    counter = 0
    ilat= BlockTools.GetIndex(self.main_dataset,"lat",str(point_coords[0]))
    ilon= BlockTools.GetIndex(self.main_dataset,"lon",str(point_coords[1]))
    #print(ilat,ilon,str(point_coords[0]),str(point_coords[1]))
    #computing composite
    for t in np.arange(0,len(self.main_dataset.time.values)):
      if dat2[t,ilat,ilon] == 0:
        dat1[t,:,:] = 0
        if additional != "":
          dat3[t,:,:] = 0
        counter += 1
    #update dat1
    dat1 = np.ma.masked_equal(dat1,0)#mask array1 when ==0
    dat1 = dat1.mean(0)
    #Compute percentage
    perc = (1-(counter/self.main_dataset.time.values.shape[0]))*100


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
      dat3 = np.ma.masked_equal(dat3,0)
      dat3 = dat3.mean(0) - meanT #temp anomaly
      max = np.amax(np.absolute(dat3))
      print(max)
      c3_range = np.linspace(-max,+max,20,endpoint = True)
      c3 = ax.contourf(lons, lats,dat3,c3_range,cmap="coolwarm", transform=datacrs) #first ploot (float)
      cbar = plt.colorbar(c3, orientation='horizontal', pad=0, label = colorbar_label ,ticks = c3_range[1:-1], aspect=40)
      cbar.ax.tick_params(labelsize=7)

    # Make some nice titles for the plot (one right, one left)
    plt.title(self.plot_title, loc='left')
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
  def PlotTracking(self,
           output,
           starting_day="2018-01-01T09:00:00.000000000",
           pIB_tracked="pIB_tracked", #name of IB matrix variable, may vary from file to file
           zg="zg", #name of zg500 matrix, may vary from file to file
           mapcrs = ccrs.PlateCarree(),
           datacrs = ccrs.PlateCarree(),
           extent =[-180,180,0,90]):

    lons,lats,dat1,dat2 = self.check_zg_pIB("zg",pIB_tracked)

    #selecting data from the desired days
    index1 = BlockTools.GetIndex(self.main_dataset,"time",starting_day)
    index2 = index1 + 15
    dat1 = dat1[index1:index2,:,:]
    dat2 = dat2[index1:index2,:,:]
    #multiple plots for having both contours and flags over different days
    dat2 = BlockTools.OrderIndexes(dat2)
    min = np.amin(dat2) #store min and max for following calculation
    max = np.amax(dat2)
    dat2 = np.ma.masked_equal(dat2,0) #mask array1 when ==0
    for i in range(0,15):
      ax = plt.subplot(5,3,i+1, projection=mapcrs, aspect = 1.0)
      ax.set_extent(extent, datacrs)
      ax.coastlines()
      #define ranges
      cs_range = np.arange(4500,6500,80) #usually it works
      #plot contour
      #dat2 is normalized in 0,1
      cs_b = ax.contour(lons, lats,dat1[i,:,:], cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#sec>    plt.clabel(cs_b, colors ="black", fontsize = 5, inline ="false")
      cb = ax.pcolor(lons,lats,dat2[i,:,:],cmap="tab20",vmin = 1, vmax = max+1,transform=datacrs)#blocking plot
      # Make some nice titles for the plot
      plt.title(starting_day[0:8] + starting_day[8:10] + ' + ' + str(i) + ' days',fontsize = 12, loc='left')

    #Export image
    try:
      plt.savefig(output,bbox_inches='tight',dpi=250)
      return 0
    except:
      print("Error code 1: no output file was given or something went wrong with matplotlib")
      return 1
