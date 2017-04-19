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

def make_plot(filename1, filename2, filename3, filename4, count):
    print filename1, filename2, filename3, filename4
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
    triang = tri.Triangulation(x, y, triangles=triang) 
    static.close()

    data1 = Dataset(filename1)
    data2 = Dataset(filename2)
    data3 = Dataset(filename3)
    data4 = Dataset(filename4)

    refl1 = np.compress(mask, data1.variables["refl10cm_1km"][0]) 
    refl2 = np.compress(mask, data2.variables["refl10cm_1km"][0])
    refl3 = np.compress(mask, data3.variables["refl10cm_1km"][0])
    refl4 = np.compress(mask, data4.variables["refl10cm_1km"][0])
    
    time = data1.variables["xtime"][0]
    time = ''.join(time)
    time = time.strip(" ")
    
    data1.close()
    data2.close()
    data3.close()
    data4.close()
    
    m.drawcoastlines(color="w")
    m.drawcountries(color="w")
    m.drawstates(color="w")
    m.drawcounties(color="w")
    m.drawmapboundary(fill_color="k")

    plot_shapefiles(m)

    cm1 = ax.tricontourf(triang, refl3, levels=np.arange(35,45,5), cmap=plt.get_cmap("Purples"), extend="max", label="Day 3")
    cm2 = ax.tricontourf(triang, refl2, levels=np.arange(35,45,5), cmap=plt.get_cmap("Greens"), extend="max", label="Day 2")
    cm3 = ax.tricontourf(triang, refl1, levels=np.arange(35,45,5), cmap=plt.get_cmap("Reds"), extend="max", label="Day 1")
    cm4 = ax.tricontourf(triang, refl4, levels=np.arange(35,45,5), cmap=plt.get_cmap("Blues"), extend="max", label="Day 4")

    plt.legend()
    
    plt.tight_layout()

    plt.savefig("/glade/u/home/khalbert/scripts/OUN_paint" + str(count).zfill(3) + ".png", pad_inches=0.1, bbox_inches="tight")
    plt.clf()
    plt.close(fig)
    end = datetime.datetime.utcnow()
    print end - start

if __name__ == "__main__":
    
    dirs = sorted(glob.glob("/glade/scratch/mpasrt/conv/2016*")) 
    print dirs
    datas1 = sorted(glob.glob(dirs[-4] + "/diag*"))
    datas2 = sorted(glob.glob(dirs[-3] + "/diag*"))
    datas3 = sorted(glob.glob(dirs[-2] + "/diag*"))
    datas4 = sorted(glob.glob(dirs[-1] + "/diag*"))

    pool = Pool(processes=25)
    counts = np.arange(0, 25, 1)
    print datas1[72 + 24]
    for count in counts[::-1]: 
        pool.apply_async(make_plot, [datas4[count], datas3[count + 24], datas2[count + 48], datas1[count + 72], count])
    #    make_plot(datas4[count], datas3[count + 24], datas2[count + 48], datas1[count + 72],  count)
    #    break
    pool.close()
    pool.join() 
