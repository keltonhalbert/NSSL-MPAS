#!/glade/u/home/khalbert/anaconda/bin/python
import numpy as np
import sharppy
from sharppy.io import pecan_decoder as dc
import glob
import datetime
from netCDF4 import Dataset
import netCDF4 as nc4

sites = {}
fcst_runs = sorted(glob.glob("/glade/p/work/khalbert/soundings/2015*"))

def find_snds():
    """
    Loop through the available soundings and store them by site.
    This is to make it easier to run through and conduct the verification
    for a given site.
    """
    for fcst in fcst_runs:
        snds = sorted(glob.glob(fcst + "/*.txt"))
        for snd in snds:
            site = snd.split("/")[-1].split(".")[0]
            try:
                sites[site].append(snd)
            except:
                sites[site] = [snd]
    return

def verify(site_id, t_idx):
    """
    Perform a forecast verification at a given site_id, date_idx, and 
    t_idx (time). Since the model is only run once daily at 00Z,
    and the model runs 5 days, 0 <= t_idx <=5.
    """
    if t_idx < 0 or t_idx > 5:
        print "Invalid verification time"
        return

    outfile = Dataset("sounding_data.nc", "a")
    ids = nc4.chartostring(outfile.variables["site_id"][:])
    site_idx = np.where(ids == site_id)[0][0]
    print site_idx
    ## convert they day x forecast to an hour forecast
    dt = 24 * t_idx
    dt = datetime.timedelta(hours=dt)
    ## select the sonding files at the site
    snds = sites[site_id]
    
    ## this stores the profile collections
    ## for each forecast day
    prof_colls = [] 
    ## loop over the soundings
    for snd in snds:
        prof_colls.append(dc.PECANDecoder(snd)._parse())
    
    ## loop over all forecast days
    ## keep track of what day we are
    cur_idx = 0
    for prof_coll in prof_colls:
        ## we're looking to verify a particular forecast time
        fcst_date = prof_coll.getCurrentDate() + dt
 
        ## loop over the current profile collection
        ## until the appropriate sounding is found
        current_date = prof_coll.getCurrentDate()
        print "Parsing Forecast Initialized: ", current_date
        print "Forecast Time to Verify: ", fcst_date
        loop_start = current_date
        while current_date != fcst_date:
            current_date = prof_coll.advanceTime(1)
            print "\tIncrementing... ", current_date
            if current_date == loop_start: break
        if current_date != fcst_date: 
            cur_idx +=1
            continue
        ## now that we're at the right date, get the profile
        fcst_prof = prof_coll.getCurrentProfs()["DETERMINISTIC"]

        ## not all forecasts have an analysis to verify against.
        ## these ensure that if they aren't found, the code doesn't
        ## crash
        found = False        
        try:
            anly_prof = prof_colls[cur_idx + t_idx].getCurrentProfs()["DETERMINISTIC"]
            found = True
        except: pass
        
        if not found:
            continue
        else:
            ## basic, non-derived variables
            tmpc_bias = fcst_prof.tmpc - anly_prof.tmpc
            dwpc_bias = fcst_prof.dwpc - anly_prof.dwpc
            ## model data is at constant height coords, so hght bias is 0
            pres_bias = fcst_prof.pres - anly_prof.pres
            uwin_bias = fcst_prof.u - anly_prof.u
            vwin_bias = fcst_prof.v - anly_prof.v

            ## derived variables
            ## pos/neg buoyancy
            mucape_bias = fcst_prof.mupcl.bplus - anly_prof.mupcl.bplus
            sfcape_bias = fcst_prof.sfcpcl.bplus - anly_prof.sfcpcl.bplus
            mucinh_bias = fcst_prof.mupcl.bminus - anly_prof.mupcl.bminus
            sfcinh_bias = fcst_prof.sfcpcl.bminus - anly_prof.sfcpcl.bminus

            ## parcel levels
            mulcl_bias = fcst_prof.mupcl.lclhght - anly_prof.mupcl.lclhght
            mulfc_bias = fcst_prof.mupcl.lfchght - anly_prof.mupcl.lfchght
            if mulfc_bias is np.ma.masked: mulfc_bias = 9.969209968386869e+36
            sflcl_bias = fcst_prof.sfcpcl.lclhght - anly_prof.sfcpcl.lclhght
            sflfc_bias = fcst_prof.sfcpcl.lfchght - anly_prof.sfcpcl.lfchght
            if sflfc_bias is np.ma.masked: sflfc_bias = 9.969209968386869e+36

            ## effective inflow layer
            eff_bot_bias = fcst_prof.ebotm - anly_prof.ebotm
            if eff_bot_bias is np.ma.masked: eff_bot_bias = 9.969209968386869e+36
            eff_top_bias = fcst_prof.etopm - anly_prof.etopm
            if eff_top_bias is np.ma.masked: eff_top_bias = 9.969209968386869e+36
            eff_depth_bias = (fcst_prof.etopm - fcst_prof.ebotm) - (anly_prof.etopm - anly_prof.ebotm)
            if eff_depth_bias is np.ma.masked: eff_depth_bias = 9.969209968386869e+36
            
            ## lapse rates & PWAT
            lr_0km_3km_bias = fcst_prof.lapserate_3km - anly_prof.lapserate_3km
            lr_3km_6km_bias = fcst_prof.lapserate_3_6km - anly_prof.lapserate_3_6km
            pwat_bias = fcst_prof.pwat - anly_prof.pwat
            
            mulfc_fcst = fcst_prof.mupcl.lfchght
            sflfc_fcst = fcst_prof.sfcpcl.lfchght
            eff_bot_fcst = fcst_prof.ebotm
            eff_top_fcst = fcst_prof.etopm
            eff_depth_fcst = fcst_prof.etopm - fcst_prof.ebotm
            mulfc_anly = anly_prof.mupcl.lfchght
            sflfc_anly = anly_prof.sfcpcl.lfchght
            eff_bot_anly = anly_prof.ebotm
            eff_top_anly = anly_prof.etopm
            eff_depth_anly = anly_prof.etopm - anly_prof.ebotm

            if mulfc_fcst is np.ma.masked: mulfc_fcst = 9.969209968386869e+36
            if sflfc_fcst is np.ma.masked: sflfc_fcst = 9.969209968386869e+36
            if eff_bot_fcst is np.ma.masked: eff_bot_fcst = 9.969209968386869e+36
            if eff_top_fcst is np.ma.masked: eff_top_fcst = 9.969209968386869e+36
            if eff_depth_fcst is np.ma.masked: eff_depth_fcst = 9.969209968386869e+36
            if mulfc_anly is np.ma.masked: mulfc_anly = 9.969209968386869e+36
            if sflfc_anly is np.ma.masked: sflfc_anly = 9.969209968386869e+36
            if eff_bot_anly is np.ma.masked: eff_bot_anly = 9.969209968386869e+36
            if eff_top_anly is np.ma.masked: eff_top_anly = 9.969209968386869e+36
            if eff_depth_anly is np.ma.masked: eff_depth_anly = 9.969209968386869e+36
            ## kinematic variables
            ## shear
            #shr_sfc_1km_bias = fcst_prof.sfc_1km_shear - anly_prof.sfc_1km_shear
            #shr_sfc_3km_bias = fcst_prof.sfc_3km_shear - anly_prof.sfc_3km_shear
            #shr_sfc_6km_bias = fcst_prof.sfc_6km_shear - anly_prof.sfc_6km_shear
            #shr_sfc_9km_bias = fcst_prof.sfc_9km_shear - anly_prof.sfc_9km_shear
            #eff_shr_bias = fcst_prof.eff_shear - anly_prof.eff_shear
            
            ## helicity
            srh_eff_bias = fcst_prof.right_esrh[0] - anly_prof.right_esrh[0]
            #srh_1km_bias = fcst_prof.right_srh1km[0] - anly_prof.right_srh1km[0]
            #srh_3km_bias = fcst_prof.right_srh3km[0] - anly_prof.right_srh3km[0]
            
            ## critical angle
            #crit_angl_bias = fcst_prof.right_critical_angle - anly_prof.right_critical_angle
            
            ## write the variables
            outfile.variables["tmpc_bias"][site_idx, :, t_idx, cur_idx] = tmpc_bias        
            outfile.variables["dwpc_bias"][site_idx, :, t_idx, cur_idx] = dwpc_bias        
            outfile.variables["uwin_bias"][site_idx, :, t_idx, cur_idx] = uwin_bias        
            outfile.variables["vwin_bias"][site_idx, :, t_idx, cur_idx] = vwin_bias
            outfile.variables["pres_bias"][site_idx, :, t_idx, cur_idx] = pres_bias
            outfile.variables["mucape_bias"][site_idx, t_idx, cur_idx] = mucape_bias
            outfile.variables["sfcape_bias"][site_idx, t_idx, cur_idx] = sfcape_bias
            outfile.variables["mucinh_bias"][site_idx, t_idx, cur_idx] = mucinh_bias
            outfile.variables["sfcinh_bias"][site_idx, t_idx, cur_idx] = sfcinh_bias
            outfile.variables["mulcl_bias"][site_idx, t_idx, cur_idx] = mulcl_bias
            outfile.variables["mulfc_bias"][site_idx, t_idx, cur_idx] = mulfc_bias
            outfile.variables["sflcl_bias"][site_idx, t_idx, cur_idx] = sflcl_bias
            outfile.variables["sflfc_bias"][site_idx, t_idx, cur_idx] = sflfc_bias; 
            outfile.variables["eff_top_bias"][site_idx, t_idx, cur_idx] = eff_top_bias
            outfile.variables["eff_bot_bias"][site_idx, t_idx, cur_idx] = eff_bot_bias
            outfile.variables["eff_depth_bias"][site_idx, t_idx, cur_idx] = eff_depth_bias
            outfile.variables["lr_0km_3km_bias"][site_idx, t_idx, cur_idx] = lr_0km_3km_bias
            outfile.variables["lr_3km_6km_bias"][site_idx, t_idx, cur_idx] = lr_3km_6km_bias
            outfile.variables["pwat_bias"][site_idx, t_idx, cur_idx] = pwat_bias

            outfile.variables["tmpc_fcst"][site_idx, :, t_idx, cur_idx] = fcst_prof.tmpc
            outfile.variables["dwpc_fcst"][site_idx, :, t_idx, cur_idx] = fcst_prof.dwpc
            outfile.variables["uwin_fcst"][site_idx, :, t_idx, cur_idx] = fcst_prof.u
            outfile.variables["vwin_fcst"][site_idx, :, t_idx, cur_idx] = fcst_prof.v
            outfile.variables["pres_fcst"][site_idx, :, t_idx, cur_idx] = fcst_prof.pres
            outfile.variables["hght_fcst"][site_idx, :, t_idx, cur_idx] = fcst_prof.hght
            outfile.variables["mucape_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.mupcl.bplus
            outfile.variables["sfcape_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.sfcpcl.bplus
            outfile.variables["mucinh_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.mupcl.bminus
            outfile.variables["sfcinh_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.sfcpcl.bminus
            outfile.variables["mulcl_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.mupcl.lclhght
            outfile.variables["mulfc_fcst"][site_idx, t_idx, cur_idx] = mulfc_fcst 
            outfile.variables["sflcl_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.sfcpcl.lclhght
            outfile.variables["sflfc_fcst"][site_idx, t_idx, cur_idx] = sflfc_fcst 
            outfile.variables["eff_top_fcst"][site_idx, t_idx, cur_idx] = eff_top_fcst
            outfile.variables["eff_bot_fcst"][site_idx, t_idx, cur_idx] = eff_bot_fcst
            outfile.variables["eff_depth_fcst"][site_idx, t_idx, cur_idx] = eff_depth_fcst
            outfile.variables["lr_0km_3km_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.lapserate_3km
            outfile.variables["lr_3km_6km_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.lapserate_3_6km
            outfile.variables["pwat_fcst"][site_idx, t_idx, cur_idx] = fcst_prof.pwat
            
            outfile.variables["tmpc_anly"][site_idx, :, t_idx, cur_idx] = anly_prof.tmpc
            outfile.variables["dwpc_anly"][site_idx, :, t_idx, cur_idx] = anly_prof.dwpc
            outfile.variables["uwin_anly"][site_idx, :, t_idx, cur_idx] = anly_prof.u
            outfile.variables["vwin_anly"][site_idx, :, t_idx, cur_idx] = anly_prof.v
            outfile.variables["pres_anly"][site_idx, :, t_idx, cur_idx] = anly_prof.pres
            outfile.variables["hght_anly"][site_idx, :, t_idx, cur_idx] = anly_prof.hght
            outfile.variables["mucape_anly"][site_idx, t_idx, cur_idx] = anly_prof.mupcl.bplus
            outfile.variables["sfcape_anly"][site_idx, t_idx, cur_idx] = anly_prof.sfcpcl.bplus
            outfile.variables["mucinh_anly"][site_idx, t_idx, cur_idx] = anly_prof.mupcl.bminus
            outfile.variables["sfcinh_anly"][site_idx, t_idx, cur_idx] = anly_prof.sfcpcl.bminus
            outfile.variables["mulcl_anly"][site_idx, t_idx, cur_idx] = anly_prof.mupcl.lclhght
            outfile.variables["mulfc_anly"][site_idx, t_idx, cur_idx] = mulfc_anly
            outfile.variables["sflcl_anly"][site_idx, t_idx, cur_idx] = anly_prof.sfcpcl.lclhght
            outfile.variables["sflfc_anly"][site_idx, t_idx, cur_idx] = sflfc_anly
            outfile.variables["eff_top_anly"][site_idx, t_idx, cur_idx] = eff_top_anly
            outfile.variables["eff_bot_anly"][site_idx, t_idx, cur_idx] = eff_bot_anly
            outfile.variables["eff_depth_anly"][site_idx, t_idx, cur_idx] = eff_depth_anly
            outfile.variables["lr_0km_3km_anly"][site_idx, t_idx, cur_idx] = anly_prof.lapserate_3km
            outfile.variables["lr_3km_6km_anly"][site_idx, t_idx, cur_idx] = anly_prof.lapserate_3_6km
            outfile.variables["pwat_anly"][site_idx, t_idx, cur_idx] = anly_prof.pwat
        ## increment the day
        cur_idx += 1
    outfile.close()
    return

def create_ncfile():
    nc = Dataset("sounding_data.nc", "w", format="NETCDF4" )
    nc.createDimension("day", len(sites["OUN"])) 
    nc.createDimension("sites", len(sites.keys()))
    nc.createDimension("forecast", 6)
    nc.createDimension("vert", 54)
    nc.createDimension("strlen", 11)

    dims_vert = ("sites", "vert", "forecast", "day")
    dims = ("sites", "forecast", "day")
    ids = nc.createVariable("site_id", "S1", ("sites", "strlen") )
    ids[:] = nc4.stringtochar(np.array(sites.keys()))

    nc.createVariable("stnLon", "f8", ("sites"))
    nc.createVariable("stnLat", "f8", ("sites"))

    nc.createVariable("tmpc_bias", "f8", dims_vert)
    nc.createVariable("dwpc_bias", "f8", dims_vert)
    nc.createVariable("uwin_bias", "f8", dims_vert)
    nc.createVariable("vwin_bias", "f8", dims_vert)
    nc.createVariable("pres_bias", "f8", dims_vert)
    nc.createVariable("mucape_bias", "f8", dims)
    nc.createVariable("sfcape_bias", "f8", dims)
    nc.createVariable("mucinh_bias", "f8", dims)
    nc.createVariable("sfcinh_bias", "f8", dims)
    nc.createVariable("mulcl_bias", "f8", dims)
    nc.createVariable("mulfc_bias", "f8", dims)
    nc.createVariable("sflcl_bias", "f8", dims)
    nc.createVariable("sflfc_bias", "f8", dims)
    nc.createVariable("eff_top_bias", "f8", dims)
    nc.createVariable("eff_bot_bias", "f8", dims)
    nc.createVariable("eff_depth_bias", "f8", dims)
    nc.createVariable("lr_0km_3km_bias", "f8", dims)
    nc.createVariable("lr_3km_6km_bias", "f8", dims)
    nc.createVariable("pwat_bias", "f8", dims)

    ## analyses and forecass
    nc.createVariable("tmpc_fcst", "f8", dims_vert)
    nc.createVariable("dwpc_fcst", "f8", dims_vert)
    nc.createVariable("uwin_fcst", "f8", dims_vert)
    nc.createVariable("vwin_fcst", "f8", dims_vert)
    nc.createVariable("pres_fcst", "f8", dims_vert)
    nc.createVariable("hght_fcst", "f8", dims_vert)
    nc.createVariable("tmpc_anly", "f8", dims_vert)
    nc.createVariable("dwpc_anly", "f8", dims_vert)
    nc.createVariable("uwin_anly", "f8", dims_vert)
    nc.createVariable("vwin_anly", "f8", dims_vert)
    nc.createVariable("pres_anly", "f8", dims_vert)
    nc.createVariable("hght_anly", "f8", dims_vert)


    nc.createVariable("mucape_fcst", "f8", dims)
    nc.createVariable("sfcape_fcst", "f8", dims)
    nc.createVariable("mucinh_fcst", "f8", dims)
    nc.createVariable("sfcinh_fcst", "f8", dims)
    nc.createVariable("mulcl_fcst", "f8", dims)
    nc.createVariable("mulfc_fcst", "f8", dims)
    nc.createVariable("sflcl_fcst", "f8", dims)
    nc.createVariable("sflfc_fcst", "f8", dims)
    nc.createVariable("mucape_anly", "f8", dims)
    nc.createVariable("sfcape_anly", "f8", dims)
    nc.createVariable("mucinh_anly", "f8", dims)
    nc.createVariable("sfcinh_anly", "f8", dims)
    nc.createVariable("mulcl_anly", "f8", dims)
    nc.createVariable("mulfc_anly", "f8", dims)
    nc.createVariable("sflcl_anly", "f8", dims)
    nc.createVariable("sflfc_anly", "f8", dims)


    nc.createVariable("eff_top_fcst", "f8", dims)
    nc.createVariable("eff_bot_fcst", "f8", dims)
    nc.createVariable("eff_depth_fcst", "f8", dims)
    nc.createVariable("lr_0km_3km_fcst", "f8", dims)
    nc.createVariable("lr_3km_6km_fcst", "f8", dims)
    nc.createVariable("pwat_fcst", "f8", dims)
    nc.createVariable("eff_top_anly", "f8", dims)
    nc.createVariable("eff_bot_anly", "f8", dims)
    nc.createVariable("eff_depth_anly", "f8", dims)
    nc.createVariable("lr_0km_3km_anly", "f8", dims)
    nc.createVariable("lr_3km_6km_anly", "f8", dims)
    nc.createVariable("pwat_anly", "f8", dims)

    nc.close()


if __name__ == "__main__":

    minlat = 34; maxlat = 37
    minlon = -100; maxlon = -95
    print "Getting sounding files..."
    find_snds()
    #create_ncfile()
    d = Dataset("sounding_data.nc")
    lons = d.variables["stnLon"][:]
    lats = d.variables["stnLat"][:]
    stid = nc4.chartostring(d.variables["site_id"][:])
    d.close()
    stn_idx = np.where(( lons <= maxlon ) & ( lons >= minlon ) & ( lats <= maxlat ) & ( lats >= minlat ))[0]
    stns = stid[:]
    print "Sounding files acquired."
    print "Conducting verification..."
    for i in np.arange(0,6):
        for stn in stns:
            print stn, i
            verify(stn, i) 
    print "Verification complete."
