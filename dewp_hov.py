from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from multiprocessing import Pool
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import glob

#MPAS data directory
MPAS_DATA_DIR = "/glade/scratch/duda/spring_exp"
STATIC = Dataset("/glade/p/work/mpasrt/rt2015/50-3/static.nc")


minlon = -105 + 360.
maxlon = -95 + 360.
minlat = 33
maxlat = 44
res = 0.025


time = np.arange(0, 6)
lons = np.arange(minlon, maxlon+res, res)
lats = np.arange(minlat, maxlat+res, res)

x, y = np.meshgrid(lons, lats)
T, X = np.meshgrid(time[::-1], lons)

cellLat = np.degrees(STATIC.variables["latCell"][:])
cellLon = np.degrees(STATIC.variables["lonCell"][:])

cellIds = np.where((cellLon >= minlon) & (cellLon <= maxlon) & (cellLat >= minlat) & (cellLat <= maxlat))[0]

runs = {}
## get all MPAS forecasts
run_dirs = sorted(glob.glob(MPAS_DATA_DIR + "/201505*"))
## loop over every forecast and place it in a dictionary
## using the day of the month (0 indexed) as the key
for i in range(len(run_dirs)):
    ## take only the 00Z runs
    runs[i] = glob.glob(run_dirs[i] + "/diagnostics*_01.00.00.nc")

print x.shape, y.shape, X.shape, T.shape

def fill_array(t_idx):
    """
    Given a forecast time step of 0-5, 0 being the 
    analysis and 5 being the day 5 forecast, compute 
    the mean bias for the month of May.
    """

    bias = []
    ## loop over every forecast run in the month of May
    for r_idx in range(len(runs)):
        try:
            ## open the forecast for that run
            fcst = Dataset(runs[r_idx][t_idx])
            print r_idx, t_idx, " FCST: " + runs[r_idx][t_idx]
        
            anl = Dataset(runs[r_idx + t_idx][0])
            print r_idx, t_idx, " ANL: " + runs[r_idx + t_idx][0]
            print "\n"
        except:
            anl = None 
        ## if there was an analysis found, compute
        if anl != None:
            dwpc_fcst = fcst.variables["dewpoint_surface"][0][cellIds]
            dwpc_anl = anl.variables["dewpoint_surface"][0][cellIds]
            print dwpc_anl.min(), dwpc_anl.max()
            
            #bias_mpas = dwpc_fcst - dwpc_anl 
            bias_mpas = dwpc_anl
            
            bias_grid = griddata((cellLon[cellIds], cellLat[cellIds]), bias_mpas, (x, y), method="linear")

            ## take the acreage across latitude 
            avg = np.nanmean(bias_grid, axis=0)
            bias.append(avg)

    bias = np.array(bias)
    mean_bias = np.nanmean(bias, axis=0)
    return mean_bias

if __name__ == "__main__":
    #pool = Pool(processes=6)
    #t_idxs = np.arange(0,6,1)
    #data = pool.map(fill_array, t_idxs)
    #dewp = np.zeros(X.shape)
    #for i in range(len(data)):
    #    dewp[:, i] = data[i]
    #np.savez("dewp_hov.npz", data=dewp)
    #"""
    datafile = np.load("dewp_hov.npz")
    data = datafile["data"][:]
    print X.shape, T.shape, data.shape
    m = Basemap(projection='lcc', resolution='i', width=1500000, height=2250000,
        lat_0=40., lon_0=-100, lat_1=30., lat_2=45.)

    plt.title("Surface Dewpoin (C)\n During HWT Spring Experiment\nAveraged from " + str(minlat) + "N to " + str(maxlat) + "N")
    plt.xlabel("Longitude")
    plt.ylabel("Forecast lead time (days)")
    xticks = np.arange(minlon, maxlon+2, 2)
    plt.xticks(xticks, xticks - 360)
    
    cf = plt.contourf(X, T, data[:, ::-1], cmap=plt.get_cmap("BrBG"), levels=np.arange(-10,31,1), extend="both")

    plt.gca().invert_yaxis()
    plt.colorbar(cf, orientation='vertical')
    plt.tight_layout()

    a = plt.axes([.65, .6, .2, .2])
    m.drawstates()
    m.plot([minlon, minlon], [minlat, maxlat], 'r-', latlon=True)
    m.plot([maxlon, maxlon], [minlat, maxlat], 'r-', latlon=True)
    m.plot([minlon, maxlon], [minlat, minlat], 'r-', latlon=True)
    m.plot([minlon, maxlon], [maxlat, maxlat], 'r-', latlon=True)
    plt.title('Region')
    plt.xticks([])
    plt.yticks([])

    plt.savefig("./sfc_dwpc_anl_hov_" + str(minlat) + "N_" + str(maxlat) + "N_mean_bias_may2015.png", bbox_inches="tight")    
    plt.close()
    #"""
