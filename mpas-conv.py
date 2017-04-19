import matplotlib as mpl
mpl.use("agg")
from mpl_toolkits.basemap import Basemap
from multiprocessing import Pool
import matplotlib.pyplot as plt
import matplotlib.tri as tri 
from netCDF4 import Dataset
import filter as mf
import numpy as np
import datetime
import glob
import os

def plot_shapefiles(m):
    states = ["KS", "NE", "CO", "OK", "SD", "AR", "TX", "MO", "MN", "IA"]
    for state in states:
        highways_info = m.readshapefile('./roads/tl_2013_'+state+'_prisecroads','highways',drawbounds=False)
        for info, highway in zip(m.highways_info, m.highways):
            #if info['RTTYP'] == 'S':
            #    x,y = zip(*highway)
            #    plt.plot(x,y, marker=None,color='#A59A2A',linewidth=0.35, zorder=2)
            if info['RTTYP'] == 'U':
                x,y = zip(*highway)
                plt.plot(x,y, marker=None,color='#A54F2A',linewidth=1.0, zorder=2)
            if info['RTTYP'] == 'I':
                x,y = zip(*highway)
                plt.plot(x,y, marker=None,color='b',linewidth=2.0, zorder=2)

def make_plot(filename, count):
    print filename 
    start = datetime.datetime.utcnow()
    fig = plt.figure(figsize=(20,20))
    ax = fig.add_subplot(111)
    #lon = -1 * (100 + count * 1.5)
    m = Basemap(projection='lcc', resolution='i', width=1500000, height=1500000,
        lat_0=40., lon_0=-100, lat_1=30., lat_2=45.) 

    static = Dataset("/glade/p/work/mpasrt/rt2016/static.nc")
    lons = np.degrees(static.variables["lonCell"][:])
    lats = np.degrees(static.variables["latCell"][:])
    x,y = m(lons, lats)
    
    cellsOnCell = static.variables["cellsOnCell"][:]
    nEdgesOnCell = static.variables["nEdgesOnCell"][:]
    nCells = 6488066
    maxEdges = 7
    fpasses = 50
    
    mask = np.logical_or(x<1.e20,y<1.e20)
    x = np.compress(mask,x)
    y = np.compress(mask,y)
 
    file = np.load("/glade/u/home/khalbert/scripts/mpas_triangles.npz")
    triang = file["plains_tri"]
    #triang2 = file["cyl_tri"]
    triang = tri.Triangulation(x, y, triangles=triang) 
    #np.savez("/glade/u/home/khalbert/scripts/mpas_triangles.npz", plains_tri=triang.triangles)
    static.close()

    data = Dataset(filename)

    cape = np.compress(mask, data.variables["refl10cm_1km"][0]) 
    #cape = np.compress(mask, data.variables["dewpoint_surface"][0])
    time = data.variables["xtime"][0]
    time = ''.join(time)
    time = time.strip(" ")
   
    u_1km = data.variables["uzonal_1km"][0]
    v_1km = data.variables["umeridional_1km"][0]
    u_6km = data.variables["uzonal_6km"][0]
    v_6km = data.variables["umeridional_6km"][0]
    
    u_shear = u_6km - u_1km
    v_shear = v_6km - v_1km
    
    shr_spd = np.sqrt((u_shear ** 2) + (v_shear ** 2)) * 1.94384449
    
    mslp2 = data.variables["cape"][0]
    mf.mpas_filter_cells(nEdgesOnCell,cellsOnCell,maxEdges,mslp2,fpasses)
    mf.mpas_filter_cells(nEdgesOnCell,cellsOnCell,maxEdges,shr_spd,fpasses)
    mslp2 = np.compress(mask, mslp2) 
    shr_spd = np.compress(mask, shr_spd)

    data.close()
    
    m.drawcoastlines(color="w")
    m.drawcountries(color="w")
    m.drawstates(color="w")
    m.drawcounties(color="w")

    m.drawmapboundary(fill_color="k")

    plot_shapefiles(m)
    cm = ax.tricontourf(triang, cape, levels=np.arange(5,75,5), cmap=plt.get_cmap("gist_ncar"), extend="max")
    cm2 = ax.tricontourf(triang, shr_spd, levels=np.arange(20,100,20), alpha=0.35, cmap=plt.get_cmap("RdPu"))
    c2 = ax.tricontour(triang, shr_spd, levels=np.arange(20,100,20), linewidths=2., linestyles="dashed", cmap=plt.get_cmap("RdPu"))
    c = ax.tricontour(triang, mslp2, levels=np.arange(500, 5500, 500), cmap=plt.get_cmap("spring_r"), linewidths=2)
    #cm = ax.tricontourf(triang, cape, levels=np.arange(4,62,2), cmap=plt.get_cmap("PuBuGn"), extend="max")
    #c = ax.tricontour(triang, mslp*.1, levels=np.arange(450, 600, 6), colors='k')
    #ax.annotate('Jamie go here!', color='r', xy=(x[max], y[max]), xytext=(x[max] + 100000, y[max] + 50000),
    #    arrowprops=dict(facecolor='red', shrink=0.05))
    
    units = "m^2/s^2"
    
    plt.clabel(c, inline=1, fmt='%1i', fontsize=20)
    plt.clabel(c2, inline=1, fmt='%1i', fontsize=20)
    plt.title("MPAS CONUS 3km Mesh\n1km AGL Reflectivity, CAPE, and 1km-6km Shear\n" + time)

    plt.tight_layout()

    plt.savefig("/glade/u/home/khalbert/scripts/OUN_refl" + str(count).zfill(3) + ".png", pad_inches=0.1, bbox_inches="tight")
    plt.clf()
    plt.close(fig)
    end = datetime.datetime.utcnow()
    print end - start

if __name__ == "__main__":
    
    dirs = sorted(glob.glob("/glade/scratch/mpasrt/conv/201605*"))
    datas = sorted(glob.glob(dirs[-2] + "/diag*"))
    pool = Pool(processes=25)
    counts = np.arange(0, len(datas), 1)

    for count in counts[::-1]: 
        pool.apply_async(make_plot, [datas[count], count])
    #    make_plot(datas[count], count)
    #    break
    pool.close()
    pool.join() 
