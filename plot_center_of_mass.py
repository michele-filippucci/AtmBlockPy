#!/usr/bin/env python
# coding: utf-8

# # Study of the dynamic of blocking
# In this notebook a plot of the blocking over europe is performed, studying the movement of the center of mass of the various events

# In[17]:


#import some basic classes
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy import stats
import scipy

#my own classes
import sys
sys.path.append("/home/guest/work/michele/prog/GitPROG/AtmBlockPy/lib")
import BlockPlots as BP
import BlockTools as BT


# Declaring the file location

# In[2]:


#coordinates of the point where we look for blocking
plon = 10.0
plat = 55.0

fn = "/home/guest/work/michele/data/ERA5/processed/"+         "ERA5_pIB_daily_djfm_northem_1979-2019.nc"
fn_out = "/home/guest/work/michele/data/ERA5/processed/"+         "ERA5_pIB_tracked_daily_djfm_northem_1979-2019.nc"
output = "/home/guest/work/michele/prog/plots/center_of_mass/traj_lag_ERA5_1979-2019_"+str(plat)+"N"+str(plon)+"E.png"


# ### Produce the array of the center of masses coordinates
# 
# In the following part the methods from BlockTools class are used.

# In[11]:


#loading dataset
ds = xr.load_dataset(fn_out)
lons,lats,dat1,dat2 = BP.check_zg_pIB(ds,"zg","pIB_tracked")
#find indexes
lonE = [BT.GetIndex(ds,"lon","-10.0"),BT.GetIndex(ds,"lon","40.0")] #Europe longitudes
latE = [BT.GetIndex(ds,"lat","30.0"),BT.GetIndex(ds,"lat","75.0")] #Europe Latitudes
#lonHotSpot = BT.GetIndex(ds,"lon",str(plon))
#latHotSpot = BT.GetIndex(ds,"lat",str(plat)) 
dat2 = np.ma.masked_equal(dat2,0)#mask array1 when ==0

"""
This loop finds the blocking events whose initial area crosses the point
(plon,plat)
"""
#create a list of indexes
indexes = []
for t in tqdm(range(len(ds["time"].values))):
  for l in np.unique(dat2[t,:,:]):
    print(l)
    if l in dat2[t,latE[0]:latE[1],lonE[0]:lonE[1]] and not l in dat2[t-1,:,:]: #box case
    #if l==dat2[t,latHotSpot,lonHotSpot] and not l in dat2[t-1,:,:]:
      indexes.append(l)
      
"""
Another loop could be implemented that ask for a condition on
the initial center of mass
"""

"""
A list of trajectories is then produced
"""

print(indexes)
xs = []
ys = []
for l in tqdm(indexes):
  x,y = BT.CenterofMass(dat2,l)
  xs.append(x)
  ys.append(y)
#print(xs,ys)


# In[12]:


# ds = xr.load_dataset(fn_out)
# starting_day="1997-01-01T09:00:00.000000000"
# lons,lats,dat1,dat2 = BP.check_zg_pIB(ds,"zg","pIB_tracked")
# #selecting data from the desired days
# index1 = BP.BlockTools.GetIndex(ds,"time",starting_day)
# index2 = index1 + 50
# tuple = dat2
# tuple = tuple[index1:index2,:,:]
# #tuple = BT.OrderIndexes(tuple)
# indexes = np.unique(tuple)
# x,y = BT.CenterofMass(tuple,2)
# #print(indexes,x,y)


# ### Trajectories Plot
# 
# In this plot only blocking event that passes through Europe are plotted.

# In[13]:


mapcrs = ccrs.PlateCarree()
datacrs = ccrs.PlateCarree()
extent =[plon-50,plon+50,plat-20,plat+20]
#extent =[-180,180,0,90]
print(indexes)


# In[14]:


ax = plt.subplot(111, projection=mapcrs)
ax.set_extent(extent, datacrs)
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.1, linestyle='--')
for i in range(len(indexes)):
  cm = plt.plot(ys[i],xs[i],transform=datacrs)
  
plt.plot(plon,plat,'or',transform=datacrs)#transorm=datacrs)
plt.title("Tracked center of mass of blocking events born in "+str(plat) +"N"+str(plon)+"E", loc = 'left')

plt.savefig(output,bbox_inches='tight',dpi=180)


# In[15]:


"""
Here we compute the distances between the last points of the trajectories
and the first points.
"""
diff = [] #vector
dist = [] #scalar distance (euclidean)
long_dist = [] #scalar horizontal distance
long_versor = np.array((0,1))
for i in range(len(indexes)):
  diff.append(np.array((xs[i][-1]-xs[i][0],ys[i][-1]-ys[i][0])))
  dist.append(np.linalg.norm(diff[i]))
  long_dist.append(np.dot(long_versor,diff[i]))
print("Mean distance travelled along longitude by a blocking event: " + str(np.mean(long_dist)))
print("Dev. st. of the above variable: " + str(np.std(long_dist)))
print("Mean distance (euclidean): " + str(np.mean(dist)))


# In[18]:


plt.hist(long_dist,bins=7,density=True)
plt.title("long dist starting point / end point blocking events")
plt.xlabel("Longitude (°E)")
plt.ylabel("Probability Density")

mean = np.mean(long_dist)
standard_deviation = np.std(long_dist) 

x_values = np.arange(-60, 60, 0.5)
y_values = scipy.stats.norm(mean, standard_deviation)
plt.plot(x_values, y_values.pdf(x_values))
plt.show


# In[ ]:




