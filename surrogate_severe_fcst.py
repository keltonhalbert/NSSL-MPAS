import glob
import numpy as np
import netCDF4 as nc4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable

m = Basemap(projection='lcc', resolution='i', width=5000000, height=2750000,
        lat_0=38., lon_0=-98, lat_1=25., lat_2=50.)

verification_grid = np.load("/glade/u/home/khalbert/scripts/CONUS_80km.npz")
lons, lats = verification_grid["lons"], verification_grid["lats"]
X, Y = m(lons, lats)

static = Dataset("/glade/p/work/mpasrt/rt2016/static.nc")
mpas_lons = static.variables["lonCell"][:]
mpas_lats = static.variables["latCell"][:]

print mpas_lons.shape, mpas_lats.shape
def create_ncfile(filename):
    d = Dataset(filename, 'w')
    
    d.createDimension("date", 31)
    d.createDimension("forecast", 6)
    d.createDimension("80km_grid_x", X.shape[1])
    d.createDimension("80km_grid_y", X.shape[0])
    d.createDimension("datelen", len("2016-05-01_00:00:00"))
    
    XLONS = d.createVariable("80km_XLONS", 'f4', ("80km_grid_y", "80km_grid_x"))
    XLATS = d.createVariable("80km_XLATS", 'f4', ("80km_grid_y", "80km_grid_x"))
    d.createVariable("dates", 'S1', ("date", "forecast", "datelen"))
    
    XLONS[:] = lons[:]
    XLATS[:] = lats[:]
    
    d.createVariable("80km_SSR_120km", 'f4', ("date", "forecast", "80km_grid_y", "80km_grid_x"))
    d.close()


def grid_all_periods(forecast_dir, date_idx):

    UH_THRESH_VAL = 100

    all_ncfiles = sorted(glob.glob(forecast_dir + "/diag*.nc"))
    period_start_end_ncfiles = sorted(glob.glob(forecast_dir + "/diag*12.00.00.nc")) 
    period_idxs = np.in1d(all_ncfiles, period_start_end_ncfiles).nonzero()[0]

    fcst_idx = 0
    for i in range(len(period_idxs) - 1):
        
        SSR = np.zeros(X.shape)
        d = Dataset("/glade/u/home/khalbert/scripts/201605_mpas_grids.nc", 'a')
        
        valid_start = period_idxs[i]
        valid_end = period_idxs[i + 1] + 1

        files = all_ncfiles[valid_start:valid_end]
        print "START: ", all_ncfiles[valid_start]
        date_data = Dataset(all_ncfiles[valid_start])
        fcst_date_str = str(nc4.chartostring(date_data.variables["xtime"][0][:])).replace(' ', '')
        fcst_date = nc4.stringtochar(np.array(fcst_date_str, dtype=str))
        print fcst_date
        date_data.close()
        print "END: ", all_ncfiles[valid_end - 1]
        print "\t FCST: ", fcst_idx

        for file in files:
            print "\t FILE: ", file
            data = Dataset(file)
            uh = data.variables["updraft_helicity_max"][0]
            mpas_storm_reports = np.where(uh >= UH_THRESH_VAL)[0]
            mpas_lon_reports = np.degrees(mpas_lons[mpas_storm_reports])
            mpas_lat_reports = np.degrees(mpas_lats[mpas_storm_reports])

            x_reports, y_reports = m(mpas_lon_reports, mpas_lat_reports)
            print "\t NUM REPORTS: ", len(x_reports)
            for x_report, y_report in zip(x_reports, y_reports):
                dist = (x_report - X)**2 + (y_report - Y)**2
                min_idx = np.where( dist == dist.min() )
                #if m.is_land(x_report, y_report): 
                SSR[min_idx[0], min_idx[1]] = 1

            data.close()
        print "\t NUM GRID POINTS: ", len( SSR[ SSR > 0 ] )
        out = gaussian_filter(SSR, sigma=1.5)
        d.variables["80km_SSR_120km"][date_idx, fcst_idx, :, :] = out[:]    
        d.variables["dates"][date_idx, fcst_idx, :] = fcst_date[:]
        d.close()
        fcst_idx += 1

if __name__ == "__main__":
    #create_ncfile("201605_mpas_grids.nc")
    dirs = sorted(glob.glob("/glade/scratch/mpasrt/conv/201605*"))
    dir = dirs[-1]
    grid_all_periods(dir, len(dirs)-1)
