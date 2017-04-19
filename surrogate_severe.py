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

#verification_grid = np.load("CONUS_80km.npz")
#lons, lats = verification_grid["lons"], verification_grid["lats"]
#X, Y = m(lons, lats)


static = Dataset("/glade/p/work/mpasrt/rt2016/static.nc")
#static = Dataset("/glade/p/work/mpasrt/rt2015/50-3/static.nc")
mpas_lons = static.variables["lonCell"][:]
mpas_lats = static.variables["latCell"][:]


def create_heat_ncfile(filename):
    d = Dataset(filename, 'w')

    lat_n = np.radians(50.)
    lat_s = np.radians(25.)
    lon_w = np.radians(-110. + 360)
    lon_e = np.radians(-85 + 360)


    idxs = np.where((mpas_lons > lon_w) & (mpas_lons < lon_e) & (mpas_lats < lat_n) & (mpas_lats > lat_s))[0]
    
    d.createDimension("date", 31)
    d.createDimension("fhours", 121)
    d.createDimension("numCells", idxs.shape[0])
    d.createDimension("datelen", len("2015-05-01_00:00:00"))
    
    #cellLon = d.createVariable("lonCell", 'f4', ("numCells"))
    #cellLat = d.createVariable("latCell", 'f4', ("numCells"))
    #d.createVariable("dates", 'S1', ("date", "fhours", "datelen"))
    
    #cellLon[:] = mpas_lons[idxs]
    #cellLat[:] = mpas_lats[idxs]

    #temp = d.createVariable("TMPC_sfc", "f4", ("date", "fhours", "numCells"))
    #dwpc = d.createVariable("DWPC_sfc", "f4", ("date", "fhours", "numCells"))
    #refl = d.createVariable("REFL_max", "f4", ("date", "fhours", "numCells"))
    #wspd = d.createVariable("WSPD_lv1", "f4", ("date", "fhours", "numCells"))   
    updr = d.createVariable("Updraft", "f4", ("date", "fhours", "numCells"))
  
    d.close()

#print mpas_lons.shape, mpas_lats.shape
def create_ncfile(filename):
    d = Dataset(filename, 'a')
    
    d.createDimension("date", 31)
    d.createDimension("forecast", 6)
    d.createDimension("80km_grid_x", X.shape[1])
    d.createDimension("80km_grid_y", X.shape[0])
    d.createDimension("datelen", len("2015-05-01_00:00:00"))
    
    XLONS = d.createVariable("80km_XLONS", 'f4', ("80km_grid_y", "80km_grid_x"))
    XLATS = d.createVariable("80km_XLATS", 'f4', ("80km_grid_y", "80km_grid_x"))
    d.createVariable("dates", 'S1', ("date", "forecast", "datelen"))
    
    XLONS[:] = lons[:]
    XLATS[:] = lats[:]
    
    for uh in np.arange(50, 575, 25):
        d.createVariable("80km_SSR_120km_UH" + str(uh), 'f4', ("date", "forecast", "80km_grid_y", "80km_grid_x"))
        d.createVariable("80km_BIN_120km_UH" + str(uh), 'f4', ("date", "forecast", "80km_grid_y", "80km_grid_x"))
    d.close()


def grid_all_periods(forecast_dir, date_idx):

    all_ncfiles = sorted(glob.glob(forecast_dir + "/diag*.nc"))
    period_start_end_ncfiles = sorted(glob.glob(forecast_dir + "/diag*12.00.00.nc")) 
    period_idxs = np.in1d(all_ncfiles, period_start_end_ncfiles).nonzero()[0]

    ## 50-575-25
    for UH_THRESH_VAL in range(250, 375, 25):
        print "UH: " + str(UH_THRESH_VAL)
        fcst_idx = 0
        for i in range(len(period_idxs) - 1):
        
            SSR = np.zeros(X.shape)
            d = Dataset("201505_mpas_grids.nc", 'a')
        
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
            d.variables["80km_SSR_120km_UH" + str(UH_THRESH_VAL)][date_idx, fcst_idx, :, :] = out[:]    
            d.variables["80km_BIN_120km_UH" + str(UH_THRESH_VAL)][date_idx, fcst_idx, :, :] = SSR[:]    
            d.variables["dates"][date_idx, fcst_idx, :] = fcst_date[:]
            d.close()
            fcst_idx += 1

def copy_fields(forecast_dir, date_idx):

    out = Dataset("201605_mpas_updraft.nc", "a")
    all_ncfiles = sorted(glob.glob(forecast_dir + "/diag*.nc"))
    period_start_end_ncfiles = sorted(glob.glob(forecast_dir + "/diag*.nc"))
    #period_idxs = np.in1d(all_ncfiles, period_start_end_ncfiles).nonzero()[0]
    #valid_start = period_idxs[0]
    #valid_end = period_idxs[1] + 1
    valid_start = 0
    valid_end = len(period_start_end_ncfiles)

    lat_n = np.radians(50.)
    lat_s = np.radians(25.)
    lon_w = np.radians(-110. + 360)
    lon_e = np.radians(-85 + 360)


    idxs = np.where((mpas_lons > lon_w) & (mpas_lons < lon_e) & (mpas_lats < lat_n) & (mpas_lats > lat_s))[0]


    count = 0
    for i in range(valid_start, valid_end):
        ncfile = all_ncfiles[i]
        print date_idx, ncfile
        in_data = Dataset(ncfile)
        
        #tmpc = in_data.variables["temperature_surface"][0]
        #print tmpc.shape, idxs.shape
        #dwpc = in_data.variables["dewpoint_surface"][0]
        #wspd = in_data.variables["wind_speed_level1_max"][0]
        #refl = in_data.variables["refl10cm_max"][0]
        updraft = in_data.variables["w_velocity_max"][0]
        
        #out.variables["TMPC_sfc"][date_idx, count, :] = tmpc[idxs]
        #out.variables["DWPC_sfc"][date_idx, count, :] = dwpc[idxs]
        #out.variables["WSPD_lv1"][date_idx, count, :] = wspd[idxs]
        #out.variables["REFL_max"][date_idx, count, :] = refl[idxs]
        out.variables["Updraft"][date_idx, count, :] = updraft[idxs]
        count += 1
    out.close()


if __name__ == "__main__":
    #create_ncfile("201605_mpas_grids.nc")
    create_heat_ncfile("201605_mpas_updraft.nc")
    #dirs = sorted(glob.glob("/glade/scratch/mwong/mpas-wrf/runs/201605*/mpas"))
    #dirs = sorted(glob.glob("/glade/scratch/duda/spring_exp/201505*"))
    dirs = sorted(glob.glob("/glade/scratch/mpasrt/conv/201605*"))
    for i in range(len(dirs)):
        dir = dirs[i]
        print dir
        copy_fields(dir, i)
    #for i in range(len(dirs)):
        #dir = dirs[i]
        #grid_all_periods(dir, i)

