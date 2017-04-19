from datetime import datetime
import numpy as np
import glob
import os

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
        print lines[i] 

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

if __name__ == "__main__":

    ## get the stations we will be parsing. We have to do it this way because
    ## there is a file for ever time step AND station. We want to parse
    ## by station only.
    mpas_csv = open(os.path.join(os.path.dirname(__file__), "mpas.csv"))
    csvlines = mpas_csv.readlines()
    for line in csvlines:
        line = line.split(",")
        files = sorted(glob.glob(os.path.join(os.path.dirname(__file__), "2015051300", line[-1] + "*.snd")))
        outfile = line[-1] + ".txt"
        for file in files:
            parse_file(file, outfile)
