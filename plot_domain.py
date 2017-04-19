import matplotlib as mpl
mpl.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import glob
import Nio
from mpl_toolkits.axes_grid1 import make_axes_locatable

WRF_DATA_DIR = "/glade/scratch/mwong/mpas-wrf/runs/201605*/wrf"
GFS_DATA_DIR = "/glade/scratch/khalbert/GFS/"

static = Dataset("/glade/scratch/mwong/mpas-wrf/runs/2016050100/wps/geo_em.d01.nc")
gstatic = Nio.open_file("/glade/scratch/khalbert/GFS/20160503/gfs_4_20160503_0000_000.grb2")
lons = static.variables["XLONG_M"][0]
lats = static.variables["XLAT_M"][0]
glon = gstatic.variables["lon_0"][:]; glat = gstatic.variables["lat_0"][:]
glons, glats = np.meshgrid(glon, glat)


gfs_anl_files = sorted(glob.glob(GFS_DATA_DIR + "/*/*.grb2"))

mean_wind_speed = 0
count = 0
#for file in gfs_anl_files:
#    print file
#    try:
#        data = Nio.open_file(file)
#        pres = data.variables["lv_ISBL0"][:] / 100.
#        pidx = np.where(pres == 250.)[0][0]
#        U = data.variables["UGRD_P0_L100_GLL0"][pidx, :, :]
#        V = data.variables["VGRD_P0_L100_GLL0"][pidx, :, :]
#        spd = np.sqrt(U**2 + V**2)
#        mean_wind_speed += spd
#        count += 1
#    except:
#        gfs_anl_files.pop(count)
#mean_wind_speed = mean_wind_speed / float(len(gfs_anl_files))

fig = plt.figure(figsize=(10,10))
ax = plt.gca()
m = Basemap(width=8000000,height=5000000,
                       rsphere=(6378137.00,6356752.3142),\
                       resolution='l',area_thresh=1000.,projection='lcc',\
                       lat_1=40,lat_2=30,lat_0=38.5,lon_0=-98.5)
m.drawstates()
m.drawcountries()
m.drawcoastlines()
m.drawparallels(np.arange(0.,81,10.), labels=[True,False,False,True])
m.drawmeridians(np.arange(10.,351.,10.), labels=[True,False,False,True])

#cm = m.contourf(glons, glats, mean_wind_speed, latlon=True, levels=np.arange(0,52,2), cmap=plt.get_cmap("BuPu"))

m.plot(lons[:, 0], lats[:, 0], 'b-', latlon=True, linewidth=2)
m.plot(lons[:, -1], lats[:, -1], 'b-', latlon=True, linewidth=2)
m.plot(lons[0, :], lats[0, :], 'b-', latlon=True, linewidth=2)
m.plot(lons[-1, :], lats[-1, :], 'b-', latlon=True, linewidth=2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.05)

#plt.colorbar(cm, cax=cax, label="m/s")
plt.tight_layout()
plt.savefig("domain.png", bbox_inches="tight")
