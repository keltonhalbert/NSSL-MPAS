#!/usr/bin/env python

from multiprocessing import Pool
from netCDF4 import Dataset
import numpy as np
import glob

staticfile = "/glade/p/work/mpasrt/rt2015/50-3/static.nc"

runs = []
dates = sorted(glob.glob("/glade/p/nmmm0031/201505*"))
for date in dates:
    runs.append( sorted(glob.glob(date + "/diagnostics.2015*_00.00.00.nc")) )


static = Dataset(staticfile)
cellLon = static.variables["lonCell"][:]
cellLat = static.variables["latCell"][:]

time = np.arange(0, 6)
outfile = Dataset("/glade/scratch/khalbert/mpas_mean_fields.nc", 'a', format="NETCDF4")

#outfile.createDimension("forecast_lead_time_days", len(time))
#outfile.createDimension("xlon", len(lons))
#outfile.createDimension("xlat", len(lats))
#outfile.createDimension("runs", len(runs))
#outfile.createDimension("nCells", len(cellLon))

## create the output variables
#write_vwind_bias = outfile.createVariable("grid_200hPa_vwind_bias", "f8", ("forecast_lead_time_days", "nCells"))
#write_uwind_bias = outfile.createVariable("grid_200hPa_uwind_bias", "f8", ("forecast_lead_time_days", "nCells"))
#write_temp_bias = outfile.createVariable("grid_700hPa_temp_bias", "f8", ("forecast_lead_time_days", "nCells"))
#write_hght_bias = outfile.createVariable("mslp_bias", "f8", ("forecast_lead_time_days", "nCells"))

#write_vwind_mean = outfile.createVariable("grid_200hPa_vwind_mean", "f8", ("forecast_lead_time_days", "nCells"))
#write_uwind_mean = outfile.createVariable("grid_200hPa_uwind_mean", "f8", ("forecast_lead_time_days", "nCells"))
#write_temp_mean = outfile.createVariable("grid_700hPa_temp_mean", "f8", ("forecast_lead_time_days", "nCells"))
#write_hght_mean = outfile.createVariable("mslp_mean", "f8", ("forecast_lead_time_days", "nCells"))
## map the output variables to the input variable type

data_map = {
    "umeridional_200hPa" : "grid_200hPa_vwind_bias",
    "uzonal_200hPa" : "grid_200hPa_uwind_bias",
#    "height_500hPa" : "grid_500hPa_hght_bias",
#    "temperature_500hPa" : "grid_500hPa_temp_bias"
#    "mslp" : "mslp_bias"
    }

#data_map = {
#    "umeridional_200hPa" : "grid_200hPa_vwind_mean",
#    "uzonal_200hPa" : "grid_200hPa_uwind_mean",
#    "height_700hPa" : "grid_700hPa_hght_mean",
#    "temperature_700hPa" : "grid_700hPa_temp_mean"
#    "mslp" : "mslp_mean"
#}

#outfile.close()

def fill_array(field):
    """
    Given a forecast time step of 0-5, 0 being the 
    analysis and 5 being the day 5 forecast, compute 
    the mean bias for the month of May.
    """
    times = np.arange(0,6,1)
    bias = []
    ## loop over every forecast run in the month of May
    for t_idx in times:
        day = []
        print "Time idx: ", t_idx
        for r_idx in range(len(runs)):
            print "Run idx: ", r_idx
            try:
                ## open the forecast for that run
                fcst = Dataset(runs[r_idx][t_idx])
                print r_idx, t_idx, " FCST: " + runs[r_idx][t_idx]
        
                anl = Dataset(runs[r_idx + t_idx][0])
                print r_idx, t_idx, " ANL: " + runs[r_idx + t_idx][0]
            except:
                anl = None 
            ## if there was an analysis found, compute
            if anl != None:
                print "Calculating fields..."
                field_fcst = fcst.variables[field][0]
                field_anl = anl.variables[field][0]
            
                bias_mpas = field_fcst - field_anl
                print "Field calculated, size: ", bias_mpas.shape
                #bias_vwin = griddata((cellLon, cellLat), bias_mpas, (x, y), method="linear")
   
                day.append(bias_mpas)
        bias.append(np.array(day))
    bias = np.array(bias)
    print bias.shape
    return bias

if __name__ == "__main__":

    #pool = Pool(processes=4)
    fields = data_map.keys()
    results = []
    for field in fields:
        results.append(fill_array(field))
#        break

#    results = np.array(pool.map(fill_array, fields))

#    print "Plotting Data..."
#    forecast = np.arange(0, 31, 1)
#    data = []
#    for day in forecast:
#        data.append(results[0][5])
#    data = np.array(data)
#    print data.shape

    print "Saving Data..."
    idx = 0
    for bias in results:
        for t_idx in np.arange(0,6,1):
            out_name = data_map[fields[idx]]
            output = bias[t_idx].mean(axis=0)
            print out_name, output.shape, outfile.variables[out_name][:].shape, np.nanmax(output)
            outfile.variables[out_name][t_idx, :] = output[:]
        idx += 1
    print "Saved!"

 #   pool.close()
 #   pool.join()

    outfile.close()

