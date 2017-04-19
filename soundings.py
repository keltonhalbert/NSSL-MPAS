import glob, os
import shutil
import tarfile
import numpy as np
from datetime import datetime

def copy_all_snd():
    """
    Copy the .tar.gz sounding files (and the untarred ones from the 31st)
    to a working directory for analysis
    """    
    ## get the list of tarfiles and loop over each file
    snd_tars = sorted(glob.glob("/glade/p/nmmm0031/201505*/snd.*.tar.gz"))
    for snd in snd_tars:
        ## split the path to get the individual path components
        dirs = snd.split("/")
        ## check if the directory for the particular date exists. If it does not
        ## then create the working directory
        if not os.path.exists("/glade/p/work/khalbert/soundings/" + dirs[4]):
            os.makedirs("/glade/p/work/khalbert/soundings/" + dirs[4])
        ## copy the tar file
        shutil.copyfile(snd, "/glade/p/work/khalbert/soundings/" + dirs[4] + "/" + dirs[-1])       
    ## May 31st doesn't have the files compressed, so copy the individual files
    snd_files = sorted(glob.glob("/glade/p/nmmm0031/2015053100/*.snd"))
    for snd in snd_files:
        dirs = snd.split("/")
        if not os.path.exists("/glade/p/work/khalbert/soundings/" + dirs[4]):
            os.makedirs("/glade/p/work/khalbert/soundings/" + dirs[4])
        shutil.copyfile(snd, "/glade/p/work/khalbert/soundings/" + dirs[4] + "/" + dirs[-1])


def untar_sounding_files():
    ## get the tar files
    tarfiles = sorted(glob.glob("/glade/p/work/khalbert/soundings/*/*.tar.gz"))
    ## loop over each file
    for tar in tarfiles:
        print tar
        dirs = tar.split("/")
        ## open and extract the tar files
        tar = tarfile.open(tar)
        tar.extractall(path=os.path.join("/",dirs[1], dirs[2], dirs[3], dirs[4], dirs[5], dirs[6]))
        tar.close()    

def convert_sounding_files():
    ## get the sounding files
    snd_files = sorted(glob.glob("/glade/p/work/khalbert/soundings/2015050100/*.snd"))
    for snd in snd_files:
        dirs = snd.split("/")
        print snd, os.path.join("/", dirs[1], dirs[2], dirs[3], dirs[4], dirs[5], dirs[6]) + "/" + dirs[7].split(".")[0] + ".txt"
        try:
            parse_file(snd, os.path.join("/", dirs[1], dirs[2], dirs[3], dirs[4], dirs[5], dirs[6])  + "/" + dirs[7].split(".")[0] + ".txt")
        except:
            print "File not found"

def parse_file(filename, outfile):
    """
    Take a .snd formatted file from GEMPAK,
    parse it, and write it out in the SHARPpy
    PECAN file format.
    """
    ## get the 2D arrays that must
    ## be passed to the parser
    pres = np.zeros((55), dtype=float)
    hght = np.zeros((55), dtype=float)
    tmpc = np.zeros((55), dtype=float)
    dwpc = np.zeros((55), dtype=float)
    wdir = np.zeros((55), dtype=float)
    wspd = np.zeros((55), dtype=float)

    ## split the file on newline characters
    ## and then split each line on whitespace 
    ## and remove the empty strings
    file = open(filename)
    lines = file.readlines()
    for i in range(0, len(lines)):
        lines[i] = filter(None, lines[i].strip("\n").split(" "))
        #print lines[i]

    ## get the meta information for the site
    time = lines[2][-1]
    stid = lines[2][2]
    slon = lines[3][5]
    slat = lines[3][2]
    selv = lines[3][-1]

    ## write the data from the file to the arrays
    count = 0
    for i in range(6, len(lines)):
        pres[count] = float(lines[i][0])
        hght[count] = float(lines[i][1])
        tmpc[count] = float(lines[i][2])
        dwpc[count] = float(lines[i][3])
        wdir[count] = float(lines[i][4])
        ## convert the speed from m/s to kts
        wspd[count] = float(lines[i][5]) * 1.94384449
        count += 1
    ## write the file
    write_file(outfile, time, stid, slon, slat, selv, pres, hght, tmpc, dwpc, wdir, wspd)

def write_file(outfile, time, stid, slon, slat, selv, pres, hght, tmpc, dwpc, wdir, wspd):
    ## given sounding profile data, write the file in the SHARPpy PECAN
    ## file format.

    ## we don't have any omega values
    #omga = np.ones(len(wdir)) * -999.0
    omga = np.zeros(len(wdir))

    ## write the header information for the data
    out_string = 'MEM = DETERMINISTIC\n'
    out_string += 'TIME = ' + time + '\n'
    out_string += 'STID = ' + stid + ' SLAT = ' + slat + ' SLON = ' + slon + " MALT = " + selv + "\n\n"
    out_string += "PRES, HGHT, TEMP, DEWP, WDIR, WSPD, OMGA\n"

    ## write out the profile data
    for i in range(len(tmpc)):
        out_string += ','.join(np.asarray([pres[i], hght[i], tmpc[i], dwpc[i], wdir[i], wspd[i], omga[i]], dtype=str)) + '\n'
    out_string += '\n\n'

    ## write the string to the output file
    file = open(outfile, "a")
    file.write(out_string)
    file.close()

def clean():
    files = glob.glob("/glade/p/work/khalbert/soundings/*/*.snd")
    for file in files:
        os.remove(file)


if __name__ == "__main__":
    #copy_all_snd()
    #untar_sounding_files()
    #convert_sounding_files()
    clean()
