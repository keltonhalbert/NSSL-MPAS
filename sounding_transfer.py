#!/usr/bin/env python

import glob
import os

if __name__ == "__main__":
    dates = sorted(glob.glob("/glade/scratch/mpasrt/conv/2016*"))
    print dates
    os.system("cd /glade/scratch/mpasrt/conv/ && tar -zcvf /glade/u/home/khalbert/" + dates[-1].split("/")[-1] + ".tar.gz " + dates[-1].split("/")[-1] + "/*.snd")
    os.system("cd /glade/u/home/khalbert/ && scp " + dates[-1].split("/")[-1] + ".tar.gz sharp:/data/soundings/mpas")
    os.system("rm /glade/u/home/khalbert/" + dates[-1].split("/")[-1] + ".tar.gz")
    os.system('ssh sharp "cd /data/soundings/mpas && tar -xzvf ' + dates[-1].split("/")[-1] + '.tar.gz && ./parse_mpas.py"')
