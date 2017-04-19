import glob, sys, os

dir_dates = sorted(glob.glob("/glade/p/nmmm0001/wrfrt/realtime_ensemble/ensf/POST/20*"))
print dir_dates
members = sorted(glob.glob(dir_dates[-1] + "/post/mem_*/"))
print members
for mem in members:
    files = sorted(glob.glob(mem + "/sound_mem*.nc"))
    os.system("scp -r " + mem + "/sound_mem* sharp:/data/soundings/ncarens/")
