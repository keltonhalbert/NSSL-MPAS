import numpy as np
import netCDF4 as nc
import mpas_laplace as ml
import glob

def laplacian_filter(cellNbr, field):
    for idx in range(cellNbr.shape[0]-1):
        try:
            field[idx] = field[idx] + 3./16. * (-4.0*field[idx] + 2./3.*(field[cellNbr[idx]].sum()))
        except: continue
    return field

##assumptions
## one file per forecast hour, 
## all cells get written to the same locations in the netcdf file

#get points
#filename = '/glade/p/work/mpasrt/rt2016/static.nc' #any file that has the lat lons in it unless they change order from file to file
filename = "/glade/scratch/mwong/mpas-wrf/runs/2016050100/wps/geo_em.d01.nc"
file = nc.Dataset(filename, 'r', format='NETCDF4')
#lon  = file.variables['lonCell'][:] #are these 0D or 2d? [0])
#lat =  file.variables['latCell'][:]
#cellNbrs =  file.variables['cellsOnCell'][:]
#nEdgesOnCell = file.variables['nEdgesOnCell'][:]
#cellsOnCell = file.variables['cellsOnCell'][:]
#maxEdges = 7 
lon = file.variables["XLONG_M"][0]
lat = file.variables["XLAT_M"][0]
file.close()

#mask = np.radians([360-120.,360-60.,25.,55.]) #just want the 3km cells or as close as possible to 3km
#a = np.logical_and(lon >= mask[0], lon <= mask[1])
#b = np.logical_and(lat >= mask[2], lat <= mask[3])
#c = np.where(np.logical_and(a,b))
#print(np.size(c)) #does this seem large enough?
#iadsfasdfi
#constants & arrays for storing
n = 5
nmax = 150 # max value assumed for UH magnitude
bins = np.zeros((31,121,nmax/n))
meta = np.zeros((31,121))


#days = sorted(glob.glob("/glade/scratch/duda/spring_exp/201605*"))
days = sorted(glob.glob("/glade/scratch/mwong/mpas-wrf/runs/201605*/wrf"))
# loop over the forecast days
h = 0
for day in days:
    k = 0
    forecast = sorted(glob.glob(day + "/diag*.nc"))
    for fcst in forecast: 
        #read the netcdf file
        filename = fcst  #add me
        print filename
        file = nc.Dataset(filename, 'r', format='NETCDF4')
        print file.variables.keys()
        alkdsjfha
        #variable = 'updraft_helicity_max' #change me if needed
        variable = "w_850hPa"
        #variable = 'w_velocity_max'
        #variable = 'UP_HELI_MAX' 
        uh = file.variables[variable][0]
        #print "FILTERING, ", uh.shape
        #uh = ml.mpas_filter_cells(nEdgesOnCell,cellsOnCell,maxEdges,uh,1) 
        #print "FILTER COMPLETE, ", uh.max()
        file.close()
        #bins[h,k,:] = np.histogram(uh[c],bins=range(0,nmax+n,n))[0]
        bins[h,k,:] = np.histogram(uh,bins=range(0,nmax+n,n))[0]
        if np.max(uh) > 0:
            meta[h,k] = 1. #values should be present, file not missing
        k += 1
    h += 1
#write data to a .npz file
wf = 'wrf_w850.npz'            
np.savez_compressed(''+wf, uh=bins,meta=meta)

