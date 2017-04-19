import matplotlib as mpl
mpl.use("agg")
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata
from multiprocessing import Pool
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import glob
import Nio
#import pygrib
from scipy.spatial import cKDTree

#MPAS data directory
MPAS_DATA_DIR = "/glade/scratch/duda/spring_exp/"
WRF_DATA_DIR = "/glade/scratch/mwong/mpas-wrf/runs/201605*/wrf"
GFS_DATA_DIR = "/glade/scratch/khalbert/GFS/"

STATIC = Dataset("/glade/p/work/mpasrt/rt2016/static.nc")
static = Dataset("/glade/scratch/mwong/mpas-wrf/runs/2016050100/wps/geo_em.d01.nc")
gstatic = Nio.open_file("/glade/scratch/khalbert/GFS/20160503/gfs_4_20160503_0000_000.grb2")
test = Nio.open_file("/glade/scratch/khalbert/GFS/20160511/gfs_4_20160511_1200_000.grb2")

gridlons = static.variables["XLONG_M"][0]
gridlats = static.variables["XLAT_M"][0]

#-134.028 -59.9719 19.8336 54.5265
#minlon = -134.028 #+ 360
minlon = -120
#minlon = xlon.min() 
#maxlon = -59.9719 #+ 360
maxlon = -70
#maxlon = xlon.max() 
#minlat = 19.8336
minlat = 30.0
#maxlat = 54.5265
#maxlat = 60.0
maxlat = 50.0
print maxlat
res = 0.1
#print xlon.shape
cellLat = np.degrees(STATIC.variables["latCell"][:])
cellLon = np.degrees(STATIC.variables["lonCell"][:])
xlon = static.variables["XLONG_M"][0]; xlat = static.variables["XLAT_M"][0]
glon = gstatic.variables["lon_0"][:]; glat = gstatic.variables["lat_0"][:]
glons, glats = np.meshgrid(glon, glat)

mpas_tree = cKDTree(zip(np.radians(cellLon), np.radians(cellLat)))
wrf_tree = cKDTree(zip(np.radians(xlon.flatten()), np.radians(xlat.flatten())))
gfs_tree = cKDTree(zip(np.radians(glons.flatten()), np.radians(glats.flatten())))


mpas_runs = {}
wrf_runs = {}
gfs_anls = {}
## get all MPAS forecasts
mpas_run_dirs = sorted(glob.glob(MPAS_DATA_DIR + "/201605*"))
wrf_run_dirs = sorted(glob.glob(WRF_DATA_DIR))
gfs_run_dirs = sorted(glob.glob(GFS_DATA_DIR + "/*"))
## loop over every forecast and place it in a dictionary
## using the day of the month (0 indexed) as the key
endings = ["00.00.00.nc", "06.00.00.nc", "12.00.00.nc", "18.00.00.nc",
           "00:00:00.nc", "06:00:00.nc", "12:00:00.nc", "18:00:00.nc"]
for i in range(len(mpas_run_dirs)):
    mpas_runs[i] = sorted(glob.glob(mpas_run_dirs[i] + "/diag*.nc"))
    wrf_runs[i] = sorted(glob.glob(wrf_run_dirs[i] + "/diags_*.nc"))

for i in range(len(mpas_runs)):
    for j in range(len(mpas_runs[i]))[::-1]:
        fname = mpas_runs[i][j]
        if fname[-1*len(endings[0]):] not in endings:
            mpas_runs[i].pop(j)

for i in range(len(wrf_runs)):
    for j in range(len(wrf_runs[i]))[::-1]:
        fname = wrf_runs[i][j]
        if fname[-1*len(endings[0]):] not in endings:
            wrf_runs[i].pop(j)

for i in range(len(gfs_run_dirs)):
    gfs_anls[i] = sorted(glob.glob(gfs_run_dirs[i] + "/*.grb2"))


def make_hov(mpas_files, wrf_files,  mpas_inds, wrf_inds, gfs_inds, tidx):

    hov_data = []
    count = 0
    gfs_files = []
    for idx in range(tidx, tidx + 5):
        gfiles = gfs_anls[idx]
        for gfile in gfiles:
            gfs_files.append(gfile)
    gfs_files.append(gfs_anls[idx+1][0])
    
    for i in range(len(mpas_files)):
        mpas_file = mpas_files[i]
        wrf_file = wrf_files[i]
        gfs_file = gfs_files[i]

        print wrf_file, mpas_file, gfs_file
        try:
            mpas_data = Dataset(mpas_file)
            wrf_data = Dataset(wrf_file)
            gfs_data = Nio.open_file(gfs_file)

            ## [[ 1000.   925.   850.   700.   600.   500.   400.   300.   250.   200.
            ##  150.   100.]]
            print wrf_data.variables["P_PL"][:]
            mpas_pv_theta = mpas_data.variables["umeridional_250hPa"][0]
            wrf_pv_theta = wrf_data.variables["V_PL"][0, -4, :, :]
            pres = gfs_data.variables["lv_ISBL0"][:] / 100.
            pidx = np.where(pres == 250.)[0][0]
            gfs_pv_theta = gfs_data.variables["VGRD_P0_L100_GLL0"][pidx, :, :]
        except:
            mpas_pv_theta = mpas_data.variables["umeridional_250hPa"][0]
            wrf_pv_theta = wrf_data.variables["V_PL"][0, -4, :, :]
            gfs_pv_theta = np.zeros(glons.shape) * np.nan 

        ## for v wind
        wrf_pv_theta = wrf_pv_theta * static.variables["COSALPHA"][0] + wrf_data.variables["U_PL"][0, -10, :, :] * static.variables["SINALPHA"][0]
	mpas_pv_grid = mpas_pv_theta.flatten()[mpas_inds].reshape(x.shape)
        wrf_pv_grid = wrf_pv_theta.flatten()[wrf_inds].reshape(x.shape)
        gfs_pv_grid = gfs_pv_theta.flatten()[gfs_inds].reshape(x.shape)

        pv_grid = wrf_pv_grid - gfs_pv_grid
        avg = pv_grid
        #pv_grid = gfs_pv_grid
	#avg = np.nanmean(pv_grid, axis=0) 
        hov_data.append(avg)
        count += 1
    return np.array(hov_data)

if __name__ == "__main__":

    time = np.arange(0, len(mpas_runs[0])) * 6
    lons = np.arange(minlon, maxlon+res, res)
    lats = np.arange(minlat, maxlat+res, res)

    x, y = np.meshgrid(lons, lats)
    print lons.min(), lons.max()
    T, X = np.meshgrid(time[::-1], lons)
    d, mpas_inds = mpas_tree.query(zip(np.radians(x.flatten() + 360.), np.radians(y.flatten())), k=1)
    d, wrf_inds = wrf_tree.query(zip(np.radians(x.flatten()), np.radians(y.flatten())), k=1)
    print np.max(d), np.degrees(np.max(d)), np.degrees(np.max(d)) * 111
    d, gfs_inds = gfs_tree.query(zip(np.radians(x.flatten() + 360.), np.radians(y.flatten())), k=1)
    accum_data = 0
    print maxlat
    for i in range(2, len(mpas_runs)):
        data = make_hov(mpas_runs[i], wrf_runs[i],  mpas_inds, wrf_inds, gfs_inds, i)
        data = np.swapaxes(data, 0, 1)
        accum_data += data
        print np.nanmax(data), np.nanmin(data)
    print accum_data.shape
    accum_data = np.nanmean(accum_data, axis=1)
   
    
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

    cm = m.contourf(x, y, accum_data, latlon=True, levels=np.arange(-50,55,5), cmap=plt.get_cmap("RdBu_r"), extend="both")

    m.plot(x[:, 0], y[:, 0], 'b-', latlon=True, linewidth=2)
    m.plot(x[:, -1], y[:, -1], 'b-', latlon=True, linewidth=2)
    m.plot(x[0, :], y[0, :], 'b-', latlon=True, linewidth=2)
    m.plot(x[-1, :], y[-1, :], 'b-', latlon=True, linewidth=2)


    m.plot(gridlons[:, 0], gridlats[:, 0], 'r-', latlon=True, linewidth=2)
    m.plot(gridlons[:, -1], gridlats[:, -1], 'r-', latlon=True, linewidth=2)
    m.plot(gridlons[0, :], gridlats[0, :], 'r-', latlon=True, linewidth=2)
    m.plot(gridlons[-1, :], gridlats[-1, :], 'r-', latlon=True, linewidth=2)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)

    plt.colorbar(cm, cax=cax, label="m/s")
    plt.tight_layout()
    plt.savefig("250hPa_domain_bias_wrf.png", bbox_inches="tight")
    """
    fig = plt.figure(figsize=(10,8))
    plt.grid(True)

    plt.title("250hPa 30 Day Average Meridional Wind Bias\n WRF - GFS INIT: 2016-05-" + str(i+1).zfill(2) + "_00\nAveraged from " + str(minlat) + "N to " + str(maxlat) + "N")
    plt.xlabel("Longitude")
    plt.ylabel("Forecast (hours)")
    plt.xlim(240-360, 290-360)
    xticks = np.arange(230-360, 310-360, 10)
    yticks = np.arange(0,121,6)
    #plt.xticks(xticks, xticks - 360)
    plt.yticks(yticks)
    cf = plt.contourf(X, T, accum_data[:, ::-1], cmap=plt.get_cmap("RdBu_r"), levels=np.arange(-6,6,1), extend="both")
    #plt.contour(X, T, data[:, ::-1], levels=[-10,-20], colors=['b', 'b'], linestyle='--')
    #    plt.contour(X, T, accum_data[:, ::-1], levels=[10,20], colors=['r', 'r'], linestyle='-')
    #plt.axvline(x=-60+360, color="k"); plt.axvline(x=-125 + 360, color="k")

    plt.gca().invert_yaxis()
    plt.colorbar(cf, orientation='vertical')
    plt.tight_layout()

    plt.savefig("./wrf_avg_bias_ckd_250hPa_umeridional_hov_0p1_" + str(i+1).zfill(2) + "_" + str(minlat) + "N_" + str(maxlat) + "N_corrected.png", bbox_inches="tight")    
    plt.close()
    """
