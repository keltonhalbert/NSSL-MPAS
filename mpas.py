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

def make_plot(filename, count):
    print filename 
    start = datetime.datetime.utcnow()
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    #lon = -1 * (100 + count * 1.5)
    m = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='i')

    static = Dataset("/glade/scratch/mpasrt/mpas_conv/2015051300/static.nc")
    lons = np.degrees(static.variables["lonCell"][:])
    lats = np.degrees(static.variables["latCell"][:])
    x,y = m(lons, lats)
    
    mask = np.logical_or(x<1.e20,y<1.e20)
    x = np.compress(mask,x)
    y = np.compress(mask,y)
 
    file = np.load("mpas_triangles.npz")
    triang = file["ortho_tri"]
    triang = tri.Triangulation(x, y, triangles=triang)
    static.close()

    data = Dataset(filename)

    #cape = np.compress(mask, data.variables["refl10cm_1km"][0]) 

    cape = np.compress(mask, data.variables["dewpoint_surface"][0])
    mslp = np.compress(mask, data.variables["mslp"][0]) 
    time = data.variables["xtime"][0]
    time = ''.join(time)
    time = time.strip(" ")
    #cape = np.compress(mask, np.sqrt(data.variables["umeridional_500hPa"][0]**2 + data.variables["uzonal_500hPa"][0]**2))
    #mslp = np.compress(mask, data.variables["height_500hPa"][0]) 

    data.close()
    m.drawcountries()
    m.drawcoastlines()
    m.drawstates()

    #cm = ax.tricontourf(triang, cape, levels=np.arange(5,75,5), cmap=plt.get_cmap("gist_ncar"), extend="max")
    
    cm = ax.tricontourf(triang, cape, levels=np.arange(10,35,2), cmap=plt.get_cmap("gist_earth_r"), extend="max")
    c = ax.tricontour(triang, mslp*.01, levels=np.arange(850, 1050, 4), colors='k')

    #cm = ax.tricontourf(triang, cape, levels=np.arange(4,62,2), cmap=plt.get_cmap("PuBuGn"), extend="max")
    #c = ax.tricontour(triang, mslp*.1, levels=np.arange(450, 600, 6), colors='k')

    #units = "dbZ"
    units = "C"
    #units = "m/s"

    plt.colorbar(cm, ax=ax, shrink=0.75, label=units, pad=0.02, aspect=25)
    plt.clabel(c, inline=1, fmt='%1i', fontsize=10)
    
    #plt.title("MPAS CONUS 3km Mesh\n1km AGL Reflectivity\n" + time)
    plt.title("MPAS CONUS 3km Mesh\nSurface Dewpoint and Pressure\n" + time)
    #plt.title("MPAS CONUS 3km Mesh\n500hPa Height and Wind\n" + time)

    plt.tight_layout()

    plt.savefig("test" + str(count).zfill(3) + ".png", pad_inches=0.1, bbox_inches="tight")
    plt.clf()
    plt.close(fig)
    end = datetime.datetime.utcnow()
    print end - start

if __name__ == "__main__":
    
    datas = sorted(glob.glob("/glade/scratch/mpasrt/mpas_conv/2015051400/diagnostic*"))
    pool = Pool(processes=20)
    counts = np.arange(0, 121, 1)

    for count in counts: 
        pool.apply_async(make_plot, [datas[count], count])
    #    make_plot(datas[count], count)
    #    break
    pool.close()
    pool.join() 
