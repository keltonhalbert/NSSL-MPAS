import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import matplotlib.tri as tri 
import filter as mf

m = Basemap(projection='lcc', resolution='i', width=1500000, height=1500000,
        lat_0=40., lon_0=-100, lat_1=30., lat_2=45.) 

m.drawstates()

static = Dataset("/glade/p/work/mpasrt/rt2015/static.nc")
lons = np.degrees(static.variables["lonCell"][:])
lats = np.degrees(static.variables["latCell"][:])
x,y = m(lons, lats)

cellsOnCell = static.variables["cellsOnCell"][:]
nEdgesOnCell = static.variables["nEdgesOnCell"][:]
nCells = 6488066
maxEdges = 7
fpasses = 50

dirs = sorted(glob.glob("/glade/scratch/mpasrt/mpas_conv/2015*")) 
datas = sorted(glob.glob(dirs[-2] + "/diagnostic*"))
data = Dataset(datas[-1])

field = data.variables["cape"][0][:]

print field.max(), field.mean()
mf.mpas_filter_cells(nEdgesOnCell,cellsOnCell,maxEdges,field,fpasses)

print field.max(), field.mean()

mask = np.logical_or(x<1.e20,y<1.e20)
x = np.compress(mask,x)
y = np.compress(mask,y)


file = np.load("mpas_triangles.npz")
triang = file["pecan_tri"]


cape = np.compress(mask, field) 

#x_plt = x[triang].T
#y_plt = y[triang].T

triang = tri.Triangulation(x, y, triangles=triang) 

#test = filters.gaussian_filter(cape, 3)

#triang = tri.Triangulation(x_plt, y_plt)

plt.tricontour(triang, cape, cmap=plt.get_cmap("spring_r"), levels=np.arange(500,5000,500))
plt.show()
plt.close()

