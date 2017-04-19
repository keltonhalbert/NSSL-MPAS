import matplotlib as mpl
mpl.use("agg")
from mpl_toolkits.basemap import Basemap
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.tri as tri 
from netCDF4 import Dataset
import numpy as np
import datetime
import glob
import os

def make_plot(filename, count):
    print filename 
    start = datetime.datetime.utcnow()
    fig = plt.figure(figsize=(27.30666,13.65333), dpi=75, frameon=False) 
    ax = fig.add_subplot(111)
    plt.axis("off")
    #lon = -1 * (100 + count * 1.5)
    m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90, llcrnrlon=0,urcrnrlon=360,resolution='i')

    dirs = sorted(glob.glob("/glade/scratch/mpasrt/mpas_conv/2015*"))
    static = Dataset(dirs[-1] + "/static.nc")
    lons = np.degrees(static.variables["lonCell"][:])
    lats = np.degrees(static.variables["latCell"][:])
    x,y = m(lons, lats)
    
    mask = np.logical_or(x<1.e20,y<1.e20)
    x = np.compress(mask,x)
    y = np.compress(mask,y)
 
    file = np.load("/glade/u/home/khalbert/scripts/mpas_triangles.npz")
    triang = file["cyl_tri"]
    triang = tri.Triangulation(x, y, triangles=triang)
    #np.savez("/glade/u/home/khalbert/scripts/mpas_triangles.npz", cyl_tri=triang.triangles, ortho_tri=triang2)
    static.close()

    data = Dataset(filename)

    cape = np.compress(mask, data.variables["refl10cm_1km"][0]) 

    #cape = np.compress(mask, data.variables["dewpoint_surface"][0])
    #mslp = np.compress(mask, data.variables["mslp"][0]) 
    time = data.variables["xtime"][0]
    time = ''.join(time)
    time = time.strip(" ")
    #cape = np.compress(mask, np.sqrt(data.variables["umeridional_250hPa"][0]**2 + data.variables["uzonal_250hPa"][0]**2))
    #mslp = np.compress(mask, data.variables["height_250hPa"][0]) 

    data.close()
    m.drawcountries()
    m.drawcoastlines()
    m.drawstates()

    cm = ax.tricontourf(triang, cape, levels=np.arange(5,75,5), cmap=plt.get_cmap("gist_ncar"), extend="max")
    
    #cm = ax.tricontourf(triang, cape, levels=np.arange(10,35,2), cmap=plt.get_cmap("gist_earth_r"), extend="max")
    #c = ax.tricontour(triang, mslp*.01, levels=np.arange(850, 1050, 4), colors='k')

    #cm = ax.tricontourf(triang, cape, levels=np.arange(4,72,2), cmap=plt.get_cmap("PuBuGn"), extend="max")
    #c = ax.tricontour(triang, mslp*.1, levels=np.arange(800, 1200, 12), colors='k')

    #units = "dbZ"
    #units = "C"

    #plt.clabel(c, inline=1, fmt='%1i', fontsize=10)
    
    props = dict(boxstyle='round', facecolor='gray', alpha=0.75) 
    ax.text(0.75, 0.5, "MPAS 1km AGL Reflectivity\n" + time, transform=ax.transAxes,
             fontsize=14, verticalalignment='top', horizontalalignment="center", bbox=props,
             weight="bold")

    ax.text(0.5, 0.5, "MPAS 1km AGL Reflectivity\n" + time, transform=ax.transAxes,
             fontsize=14, verticalalignment='top', horizontalalignment="center", bbox=props,
             weight="bold")

    ax.text(0.25, 0.5, "MPAS 1km AGL Reflectivity\n" + time, transform=ax.transAxes,
             fontsize=14, verticalalignment='top', horizontalalignment="center", bbox=props,
             weight="bold")


    #plt.title("MPAS CONUS 3km Mesh\n1km AGL Reflectivity\n" + time)
    #plt.title("MPAS CONUS 3km Mesh\nSurface Dewpoint and Pressure\n" + time)

    plt.tight_layout()

    plt.savefig("/glade/u/home/khalbert/scripts/SOS-refl" + str(count).zfill(3) + ".png", pad_inches=0, bbox_inches="tight", dpi=75)
    plt.clf()
    plt.close(fig)
    end = datetime.datetime.utcnow()
    print end - start

if __name__ == "__main__":
    dirs = sorted(glob.glob("/glade/scratch/mpasrt/mpas_conv/2015*")) 
    datas = sorted(glob.glob(dirs[-1] + "/diagnostic*"))
    pool = Pool(processes=20)
    counts = np.arange(0, len(datas), 1)

    for count in counts: 
        pool.apply_async(make_plot, [datas[count], count])
    #    make_plot(datas[count], count)
    #    break
    pool.close()
    pool.join()

    #os.system("convert " + os.path.dirname(__file__) + "/250_wind* " + os.path.dirname(__file__) + "/mpas-250hPa.gif")
    #os.system("scp " + os.path.dirname(__file__) + "/mpas-250hPa.gif arctic:/var/www/html")

