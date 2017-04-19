#!/usr/bin/env python
import matplotlib as mpl; mpl.use('agg')
from scipy.interpolate import griddata
from multiprocessing import Pool
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import glob

minlon = 0
maxlon = 360
minlat = 60
maxlat = 90
res = 0.25


time = np.arange(0, 6)
lons = np.arange(minlon, maxlon+res, res)
lats = np.arange(minlat, maxlat+res, res)

x, y = np.meshgrid(lons, lats)
T, X = np.meshgrid(time[::-1], lons)

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
        except:
            anl = None 
        ## if there was an analysis found, compute
        if anl != None:
            vwind_fcst = fcst.variables["umeridional_250hPa"][0]
            vwind_anl = anl.variables["umeridional_250hPa"][0]
            
            bias_mpas = vwind_fcst - vwind_anl 
            
            bias_vwin = griddata((cellLon, cellLat), bias_mpas, (x, y), method="linear")
   
            #avg = bias_vwin.mean(axis=0)
            bias.append(avg)

    bias = np.array(bias)
    print bias.shape
    #data = outfile.variables["250hPa_vwind_bias"][:bias.shape[0], t_idx, :] = bias[:]
    #mean_bias = bias.mean(axis=0)
    return bias


if __name__ == "__main__":

    data = Dataset("processed_60N_90N.nc")
    wind_bias = data.variables["250hPa_vwind_bias"][:]

    mean_bias = wind_bias.mean(axis=0)
    mean_bias = np.swapaxes(mean_bias, 0, 1)
    print mean_bias.max()
    days = wind_bias.shape[0]
    for day in range(days):

    plt.title("250hPa V Wind Bias\n During HWT Spring Experiment\nAveraged from " + str(minlat) + "N to " + str(maxlat) + "N")
    plt.xlabel("Longitude")
    plt.ylabel("Forecast lead time (days)")

    cf = plt.contourf(X, T, mean_bias[:, ::-1], cmap=plt.get_cmap("RdBu_r"), levels=np.arange(-2,2.25,0.25), extend="both")

    bias = np.swapaxes(wind_bias[day, :, :], 0, 1)
    print bias.shape
    cs = plt.contour(X, T, bias[:, ::-1], '--', levels=[-4, 0, 4], cmap=plt.get_cmap("PiYG"), linewidth=3)
    plt.clabel(cs, inline=1, fontsize=8)

    plt.gca().invert_yaxis()
    plt.colorbar(cf, orientation='vertical')
    plt.tight_layout()
    plt.savefig("./days/250hPa_" + str(day).zfill(2) + "_hov_60N_90N_0-360_mean_bias_may2015.png", bbox_inches="tight")    
    plt.close()
