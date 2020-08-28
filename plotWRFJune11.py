from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim, smooth2d)

# Open the NetCDF file
iday=12
prefix="../PMMSims/"
prefix=""
ncfile = Dataset(prefix+"wrfout_d04_2014-06-%2.2i_20:00:00"%iday)
ncfile2 = Dataset(prefix+"wrfout_d04_2014-06-%2.2i_23:48:00"%iday)
rain=ncfile2['RAINNC'][0,:,:]-ncfile['RAINNC'][0,:,:]
smooth_rain = smooth2d(rain, 5, cenweight=4)
print(smooth_rain.max())
eeLL=[35.372831, -82.369989670]
mpLL=[35.425783, -82.7573]
mvLL=[35.5198, -83.094764]
pcLL=[35.292858,-82.1707]
# Extract the pressure, geopotential height, and wind variables
p = getvar(ncfile, "pressure")
z = getvar(ncfile, "z", units="dm")
ua = getvar(ncfile, "ua", units="kt")
va = getvar(ncfile, "va", units="kt")
wspd = getvar(ncfile, "wspd_wdir", units="kts")[0,:]

# Interpolate geopotential height, u, and v winds to 500 hPa
ht_500 = interplevel(z, p, 850)
u_500 = interplevel(ua, p, 850)
v_500 = interplevel(va, p, 850)
wspd_500 = interplevel(wspd, p, 850)

# Get the lat/lon coordinates
lats, lons = latlon_coords(ht_500)

# Get the map projection information
cart_proj = get_cartopy(ht_500)

# Create the figure
fig = plt.figure(figsize=(12,9))
ax = plt.axes(projection=cart_proj)

# Download and add the states and coastlines
states = NaturalEarthFeature(category="cultural", scale="50m",
                             facecolor="none",
                             name="admin_1_states_provinces_shp")
ax.add_feature(states, linewidth=0.5, edgecolor="black")
ax.coastlines('50m', linewidth=0.8)

# Add the 500 hPa geopotential height contours

tpw=getvar(ncfile,"pw")

levels = np.arange(10., 70., 6.)
contours = plt.contourf(to_np(lons), to_np(lats), to_np(tpw),
                      levels=levels, cmap=get_cmap("rainbow"),
                       transform=crs.PlateCarree())
#plt.scatter([eeLL[1],mpLL[1],mvLL[1]],[eeLL[0],mpLL[0],mvLL[0]],s=60,transform=crs.PlateCarree(),color='pink')
plt.scatter([eeLL[1],mpLL[1],pcLL[1]],[eeLL[0],mpLL[0],pcLL[0]],s=60,transform=crs.PlateCarree(),color='pink')
plt.scatter([eeLL[1]],[eeLL[0]],s=60,transform=crs.PlateCarree(),color='olive')
#plt.clabel(contours, inline=1, fontsize=10, fmt="%i")

# Add the wind speed contours
from numpy import array
levels = 0.25*array([25, 30, 40,  60,  80,  100,  120, 140])
rain_contours = plt.contour(to_np(lons), to_np(lats), to_np(smooth_rain),
                           levels=levels,cmap='jet',
                           transform=crs.PlateCarree())

ns=100
# Add the 500 hPa wind barbs, only plotting every 125th data point.
plt.barbs(to_np(lons[::ns,::ns]), to_np(lats[::ns,::ns]),
          to_np(ua[10,::ns, ::ns]), to_np(va[10,::ns, ::ns]),
          transform=crs.PlateCarree(), length=6)


# Set the map bounds
ax.set_xlim(cartopy_xlim(ht_500))
ax.set_ylim(cartopy_ylim(ht_500))

ax.gridlines()
plt.colorbar(contours)
plt.title("TPW,Level 10 wind barbs (kt)")
plt.savefig("June%2.2i.png"%iday)
plt.show()
