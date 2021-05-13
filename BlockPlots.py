import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

class BlockPlots(object):
    num_of_BlockPlots = 0
    plot_title = "none"
    def __init__(self,plot_title):
        BlockPlots.num_of_BlockPlots += 1

    def read_main(self,filename):
        self.main_dataset = xr.load_dataset(filename)

    """
    PlotEventZ500 function produces a plot for a specific day and time of Z500 field
    and highlight where the Z500 anomaly index detected a blocking event.
    All requested variables have defaults except from output.
    """

    def PlotEventZ500(self,
                      output,
                      day,
                      pIB="pIB",
                      zg="zg",
                      #mapcrs = ccrs.NorthPolarStereo(),
                      #datacrs = ccrs.PlateCarree(),
                      lon=[-180,+180],
                      lat=[0,90]):
        mapcrs = ccrs.NorthPolarStereo(),
        datacrs = ccrs.PlateCarree(),

        try:
            zg = self.dataset["zg"]
            pIB = self.main_dataset["pIB"]
            string = "data read successfully"
        except:
            string = "zg variable and pIB were not found.\n \
Hint: use read_main() to load data with right attributes."
            print(string)
            return 0

        #multiple plots for having both contours and heatmap
        ax = plt.subplot(111, projection=mapcrs)
        ax.set_extent(lon[0],lon[1],lat[0],lat[1], dataccrs)
        ax.coastlines()

        #define ranges
        cs_range = np.arange(4500,6500,40) #usually it works
        cb_range = range(0,1)

        #plot contour
        cs = ax.contourf(lons, lats,dat1, cs_range,cmap="jet", transform=datacrs) #first ploot (float)
        cs_b = ax.contour(lons, lats,dat1, cs_range,colors="black",linewidths = 0.1 ,transform=datacrs)#second plot for contours
        plt.clabel(cs_b, colors  ="black", fontsize ="smaller", inline ="false")
        cb = ax.contour(lons,lats, dat2, cb_range,cmap = "flag", transform=datacrs)#blocking plot

        # Make some nice titles for the plot (one right, one left)
        plt.title(self.plot_title, loc='left')
        plt.colorbar(cs, orientation='horizontal', pad=0, label = "geopotential height" , aspect=50)

        try:
            plt.savefig(output)
            return 0
        except:
            print("Error code 1: no output file was given")
            return 1
