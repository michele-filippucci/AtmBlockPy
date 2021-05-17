import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

#my own class
from BlockTools import BlockTools
import BlockTools

#plot dimensions in inches
plt.rcParams["figure.figsize"] = (10,10)
np.set_printoptions(precision=2,threshold=np.inf)

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
    def check_zg_pIB(self,pIB="",zg=""):
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
            lons = np.append(lons,180-1e-6)#slightly different for avoiding contour bug
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

        lons,lats,dat1,dat2 = self.check_zg_pIB(pIB,zg)

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
        plt.clabel(cs_b, colors  ="black", fontsize = 5, inline ="false")
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
                      extent =[-180,180,20,70]):

        lons,lats,dat1,dat2 = self.check_zg_pIB(pIB,zg)

        if additional != "":
            try:
                dat3 = self.additional_dataset[additional].loc[:,5e+04,:,:].values
                print("ok")
                dat3 = cutil.add_cyclic_point(dat3)
                meanT = dat3.mean(0)
            except:
                print("Error code 1: invalid additional dataset.")
                return 1

        counter = 0
        ilat= BlockTools.GetIndex(self.main_dataset,"lat",str(point_coords[0]))
        ilon= BlockTools.GetIndex(self.main_dataset,"lon",str(point_coords[1]))
        print(ilat,ilon,str(point_coords[0]),str(point_coords[1]))
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
        ax.set_extent([-180, 175.5, 0, 90], ccrs.PlateCarree())
        ax.coastlines()

        #plot contour
        cs_range = np.arange(4500,6500,50)
        cs_b = ax.contour(lons,lats,dat1, cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#second plot for contours
        plt.clabel(cs_b, colors  ="black", fontsize ="smaller", inline ="false")

        #update dat3 if present
        if additional != "":
            dat3 = np.ma.masked_equal(dat3,0)
            dat3 = dat3.mean(0) - meanT #temp anomaly
            cs = ax.contourf(lons, lats,dat3,cmap="coolwarm", transform=datacrs) #first ploot (float)
            plt.colorbar(cs, orientation='horizontal', pad=0, label = "Temperature (K)" , aspect=40)

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
