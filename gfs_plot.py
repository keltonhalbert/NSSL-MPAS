import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import sys, os, glob
import numpy as np
import Nio

days = sorted(glob.glob("/glade/scratch/khalbert/GFS/201505*"))

def make_plot(idx):
    
    fig = plt.figure(figsize=(20,20), dpi=100, frameon=False) 
    ax = fig.add_subplot(111)
    plt.axis("off")
    m = Basemap( projection='npstere', boundinglat=25, lon_0=-90, resolution='c' )
    static = Nio.open_file("/glade/scratch/khalbert/GFS/2015050100/gfs_4_20150501_0000_000.grb2")
    ## lon_0 - lon
    ## lat_0 - lat
    ## lv_ISBL0 - pres
    ## HGT_P0_L100_GLL0 - height
    lon = static.variables["lon_0"][:]
    lat = static.variables["lat_0"][:]
    lons, lats = np.meshgrid(lon, lat)
    print lons.shape, lats.shape, lon.shape, lat.shape
    X, Y = m(lons, lats)

    mean_bias = []
    mean_fcst = []

    for day in days:
        file1 = glob.glob(day + "/*" + str(idx*24).zfill(3) + ".grb2")
        file2 = glob.glob(day + "/*f" + str(idx*24).zfill(3))
        if len(file1) > len(file2):
            ffile = file1[0]
            fcst = Nio.open_file(ffile)
        elif len(file2) > len(file1):
            ffile = file2[0]
        try:
            file3 = glob.glob(days[days.index(day) + idx] + "/*000.grb2")
            file4 = glob.glob(days[days.index(day) + idx] + "/*f000")
        except:
            file3 = []
            file4 = []

        if len(file3) > len(file4):
             afile = file3[0]
        elif len(file4) > len(file3):
             afile = file4[0]
        else:
             afile = None

        if afile != None:
            fcst = Nio.open_file(ffile, format="grib2")
            anly = Nio.open_file(afile, format="grib2")
            p_idx = np.where(fcst.variables["lv_ISBL0"][:] == 50000.)[0][0]
            print "FCST: " + ffile
            print "ANLY: " + afile

            day_hght_fcst = fcst.variables["UGRD_P0_L100_GLL0"][p_idx, :, :]
            day_hght_anly = anly.variables["UGRD_P0_L100_GLL0"][p_idx, :, :]
            day_hght_bias = day_hght_fcst - day_hght_anly
            mean_bias.append(day_hght_bias)
            mean_fcst.append(day_hght_fcst)

    mean_bias = np.mean(mean_bias, axis=0)
    mean_fcst = np.mean(mean_fcst, axis=0)
    mean_anl = mean_fcst - mean_bias

    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color="#727678", alpha=0.15)
    print X.shape, Y.shape, mean_bias.shape
    cm = ax.contourf(X, Y, mean_bias, levels=np.arange(-10,11,1), cmap=plt.get_cmap("RdBu_r"), extend="both")
    c1 = ax.contour(X, Y, mean_fcst, levels=np.arange(0, 51, 10), linewidths=2.0, colors='k')
    #c2 = ax.contour(X, Y, mean_anl, levels=np.arange(0, 51, 10), linewidths=2.0, colors='k', linestyles="dashed")

    #plt.clabel(c1, inline=1, fmt='%1i', fontsize=18)
    #plt.clabel(c2, inline=1, fmt='%1i', fontsize=18)


    cb = plt.colorbar(cm, orientation="vertical", pad=0, shrink=0.85)
    cb.ax.tick_params(labelsize=18)
    cb.set_label("m/s")

    labl_anly = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--', linewidth=2)
    labl_fcst = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-', linewidth=2)
    plt.legend([labl_anly, labl_fcst], ['Analysis', 'Forecast'], fontsize=24)

    props = dict(boxstyle='round', facecolor='gray', alpha=0.75)
    ax.text(0.5, 0.95, "GFS Mean Day " + str(idx) + " Forecast\n U Zonal Anl, U Zonal Fcst, U Zonal Bias\n(dm)", transform=ax.transAxes,
             fontsize=14, verticalalignment='top', horizontalalignment="center", bbox=props,
             weight="bold")

    plt.tight_layout()

    plt.savefig("/glade/u/home/khalbert/scripts/npstere-gfs-mean-uwin-500_day" + str(idx).zfill(2) + ".png", pad_inches=0, bbox_inches="tight", dpi=75)
    plt.close()
if __name__ == "__main__":
    for i in np.arange(0,6,1):
        make_plot(i)
        #break
