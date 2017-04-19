import matplotlib as mpl
mpl.use("agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.tri as tri 
from netCDF4 import Dataset
import filter as mf
import numpy as np

def make_plot(filename, t_idx):


    fig = plt.figure(figsize=(20,20), dpi=100, frameon=False) 
    ax = fig.add_subplot(111)
    plt.axis("off")
    mpl.rcParams['contour.negative_linestyle']= 'solid'
    #lon = -1 * (100 + count * 1.5)
    #m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90, llcrnrlon=0,urcrnrlon=360,resolution='i')
    m = Basemap( projection='npstere', boundinglat=25, lon_0=-90, resolution='c' )

    staticfile = "/glade/p/work/mpasrt/rt2015/50-3/static.nc"

    static = Dataset(staticfile)
    lons = np.degrees(static.variables["lonCell"][:])
    lats = np.degrees(static.variables["latCell"][:])
    x,y = m(lons, lats)
    
    cellsOnCell = static.variables["cellsOnCell"][:]
    nEdgesOnCell = static.variables["nEdgesOnCell"][:]
    nCells = 6848514
    maxEdges = 7
    fpasses = 50    
    
    def filter(field):
        return mf.mpas_filter_cells(nEdgesOnCell, cellsOnCell, maxEdges, field, fpasses)

    mask = np.logical_or(x<1.e20,y<1.e20)
    x = np.compress(mask,x)
    y = np.compress(mask,y)
 
    file = np.load("/glade/u/home/khalbert/scripts/mpas_triangles.npz")
    triang = file["npstere"]
    triang = tri.Triangulation(x, y, triangles=triang)
    #np.savez("/glade/u/home/khalbert/scripts/mpas_triangles.npz", cyl_tri=triang2, npstere=triang.triangles)
    static.close()

    data = Dataset(filename)
    #temp = np.compress(mask, data.variables["grid_500hPa_temp_mean"][t_idx]) 
    #hght = np.compress(mask, data.variables["grid_500hPa_hght_mean"][t_idx])
    #mslp = np.compress(mask, filter(data.variables["mslp_mean"][t_idx]))
    #cape = np.compress(mask, filter(data.variables["cape_mean"][t_idx]))
    uwin = np.compress(mask, data.variables["grid_200hPa_uwind_mean"][t_idx])
    #temp_bias = np.compress(mask, data.variables["grid_500hPa_temp_bias"][t_idx])
    #hght_bias = np.compress(mask, data.variables["grid_500hPa_hght_bias"][t_idx]) 
    #vwin_bias = filter(data.variables["grid_500hPa_vwind_bias"][t_idx])
    uwin_bias = data.variables["grid_200hPa_uwind_bias"][t_idx]
    #vwin_bias = np.compress(mask, vwin_bias)
    uwin_bias = np.compress(mask, uwin_bias)
    #mslp_bias = np.compress(mask, data.variables["mslp_bias"][t_idx])
    #cape_bias = np.compress(mask, data.variables["cape_bias"][t_idx])

    #hght_anl = hght - hght_bias
    uwin_anl = uwin - uwin_bias
    #print temp.max(), hght.max(), vwin.max()
   
    data.close()
    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents(color="#727678", alpha=0.15)

    print uwin.mean(), uwin.min(), uwin.max()
    ## -50, 52, 2
    #cm = ax.tricontourf(triang, hght_bias*.1, levels=np.arange(-6,6.5,0.5), cmap=plt.get_cmap("RdBu_r"), extend="both")
    #cm = ax.tricontourf(triang, temp, levels=np.arange(210,270,2), cmap=plt.get_cmap("RdYlBu_r"), extend="both")
    #cm = ax.tricontourf(triang, mslp_bias*.1, levels=np.arange(-60, 61, 1), cmap=plt.get_cmap("RdBu_r"), extend="both")
    cm = ax.tricontourf(triang, uwin_bias, levels=np.arange(-10,11,1), cmap=plt.get_cmap("RdBu_r"), extend="both")
    #c1 = ax.tricontour(triang, hght*.1, levels=np.arange(200, 600, 6), linewidths=2.0, colors='k')
    #c2 = ax.tricontour(triang, hght_anl *.1, levels=np.arange(200, 600, 6), linewidths=2.0, colors='k', linestyles="dashed")
    #c3 = ax.tricontour(triang, mslp*.01, levels=np.arange(850, 1100, 2), colors="#4b0082", linewidths=2.0)
    #c4 = ax.tricontour(triang, (mslp - mslp_bias)*.01, levels=np.arange(850, 1100, 2), colors="#4b0082", linestyles="dashed", linewidths=2.0)
    #c2 = ax.tricontour(triang, vwin_bias, levels=[-4,-2,0,2,4], cmap=plt.get_cmap("RdBu_r"), linewidths=2)
    c1 = ax.tricontour(triang, uwin, levels=np.arange(-10,51,10), linewidths=2.0, colors="k")
    c2 = ax.tricontour(triang, uwin_anl, levels=np.arange(-10,51,10), linewidths=2.0, colors="k", linestyles="dashed")

    plt.clabel(c1, inline=1, fmt='%1i', fontsize=18)
    plt.clabel(c2, inline=1, fmt="%1i", fontsize=18)
    #plt.clabel(c3, inline=1, fmt="%1i", fontsize=18)
    #plt.clabel(c4, inline=1, fmt="%1i", fontsize=18)
    cb = plt.colorbar(cm, orientation="vertical", pad=0, shrink=0.85)    
    cb.ax.tick_params(labelsize=18)
    cb.set_label('m/s')

    labl_anly = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--', linewidth=2)
    labl_fcst = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-', linewidth=2)
    plt.legend([labl_anly, labl_fcst], ['Analysis', 'Forecast'], fontsize=24)

    props = dict(boxstyle='round', facecolor='gray', alpha=0.75)
    ax.text(0.5, 0.95, "MPAS Mean 200hPa Day " + str(t_idx) + " Forecast\n U Zonal Anl, U Zonal Fcst, U Zonal Bias\n(m/s)", transform=ax.transAxes,
             fontsize=14, verticalalignment='top', horizontalalignment="center", bbox=props,
             weight="bold")


    plt.tight_layout()

    plt.savefig("/glade/u/home/khalbert/scripts/npstere-mpas-mean-200_uwin_day" + str(t_idx).zfill(2) + ".png", pad_inches=0, bbox_inches="tight", dpi=75)

if __name__ == "__main__":
    for i in np.arange(0,6,1):
        make_plot("/glade/scratch/khalbert/mpas_mean_fields.nc", i)
        #break
