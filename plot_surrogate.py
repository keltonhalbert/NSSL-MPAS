import matplotlib as mpl
mpl.use("agg")
import glob
import numpy as np
import netCDF4 as nc4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.path import Path
from mpl_toolkits.basemap import Basemap
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_day(filename="201605_wrf_grids.nc", day_of_month=0, fcst_idx=0, reports=False):
    """
    Plots the surrogate severe forecast.
    """
    csv_data = np.genfromtxt("./rpts/1605" + str(day_of_month+1).zfill(2) + "_rpts.csv", delimiter=",", skip_header=1)
    lons = csv_data[:,-2]
    lats = csv_data[:,-3]

    d = Dataset(filename)
    dates = d.variables["dates"][:]
    init = str(nc4.chartostring(dates[day_of_month, 0, :]))
    fcst = str(nc4.chartostring(dates[day_of_month, fcst_idx, :]))
    print init, fcst
    ## day of month, forecast day, gridpoints - order of  array
    data = d.variables["80km_SSR_120km_UH75"][day_of_month, fcst_idx, :, :]
    print data.max(), data.min()
    #print init, fcst
    sigma = 120

    plt.figure(figsize=(11, 8.5))
    ax = plt.gca()
    plt.title("INIT: " + init + "\nFCST: " + fcst + "\n MPAS Severe Probabilities")

    m.drawstates()
    m.drawcountries()
    m.drawcoastlines()
    m.drawlsmask(land_color='#C9C9C9',ocean_color='#90C3D4', alpha=0.25)

    cm = m.contourf(X, Y, data, vmin=0, vmax=1, levels=[.05, 0.15, 0.3, 0.45, 0.60, 0.75],
                cmap=plt.get_cmap("magma"), extend="max", apha=0.25)
    m.plot(lons, lats, 'co', latlon=True)        
        
    #create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.05)

    plt.colorbar(cm, cax=cax, label="Probability")
    plt.tight_layout()
    plt.savefig(init[0:11] + "day" + str(fcst_idx + 1) + "_sigma_" + str(int(sigma)) + "km_UH100.png",  bbox_inches="tight")
    plt.close()



def plot_day_comp(filename="201605_mpas_grids.nc", day_of_month=0):
    """
    Plots the surrogate severe forecast.
    """
    csv_data = np.genfromtxt("./rpts/1605" + str(day_of_month+1).zfill(2) + "_rpts.csv", delimiter=",", skip_header=1)
    lons = csv_data[:,-2]
    lats = csv_data[:,-3]
    d = Dataset(filename)
    dates = d.variables["dates"][:]
    #print init, fcst
    sigma = 120
    
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,12))
    init_d = day_of_month + 1
    for i in np.arange(0,4,1):
        
        subplot_idx = i + 1
        ax = axes.flatten()[i]
        ## day of month, forecast day, gridpoints - order of  array
        data = d.variables["80km_SSR_120km"][day_of_month, i, :, :]
        
        ## date strings
        init = str(nc4.chartostring(dates[day_of_month, 0, :]))
        fcst = str(nc4.chartostring(dates[day_of_month, i, :]))
        
        ax = plt.subplot(2, 2, subplot_idx)
        plt.title("INIT: " + init + "\nFCST: " + fcst + "\n Day " + str(i+1) + " MPAS Storm Reports")

        m.drawstates()
        m.drawcountries()
        m.drawcoastlines()
        m.drawlsmask(land_color='#C9C9C9',ocean_color='#90C3D4', alpha=0.25)

        cm = m.contourf(X, Y, data, vmin=0, vmax=1, levels=[.05, 0.15, 0.3, 0.45, 0.60, 0.75],
                cmap=plt.get_cmap("magma"), extend="max", apha=0.25)
        m.plot(lons, lats, 'co', latlon=True)        
        #create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        day_of_month -= 1
        
    #divider = make_axes_locatable(ax1)
    #cax = divider.append_axes("right", size="2%", pad=0.05)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([1, 0.15, 0.015, 0.7])
    fig.colorbar(cm, cax=cbar_ax, label="Probability")

    #plt.colorbar(cm, cax=cax, label="Probability")
    plt.tight_layout()
    plt.savefig(fcst[0:11] + "day1_comp_sigma_" + str(int(sigma)) + "km_UH100.png",  bbox_inches="tight")
    plt.close()
        


if __name__ == "__main__":
    dirs = sorted(glob.glob("/glade/scratch/mpasrt/conv/201605*"))
    m = Basemap(projection='lcc', resolution='i', width=5000000, height=2750000,
        lat_0=38., lon_0=-98, lat_1=25., lat_2=50.)

    verification_grid = np.load("CONUS_80km.npz")
    lons, lats = verification_grid["lons"], verification_grid["lats"]
    X, Y = m(lons, lats)
    day = len(dirs) - 1
    for day in range(0,31,1):
        plot_day(day_of_month=day, fcst_idx=0)
        plot_day(day_of_month=day, fcst_idx=1)
        plot_day(day_of_month=day, fcst_idx=2)
        plot_day(day_of_month=day, fcst_idx=3)
        #plot_day_comp(day_of_month=day)
   
