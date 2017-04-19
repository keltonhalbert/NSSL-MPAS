import glob, os

files = sorted(glob.glob("/glade/p/nmmm0001/wrfrt/realtime_ensemble/ensf/POST/2015051*/post/mem_*/diags_d02*"))

for file in files:
    paths = file.split("/")
    member = paths[-2]
    date = paths[-4]
    fname = paths[-1]

    print date, member, fname

    try:
        os.mkdir("/glade/scratch/khalbert/ncarens/" + date)
    except: pass

    try:
        os.mkdir("/glade/scratch/khalbert/ncarens/" + date + "/" + member)
    except: pass
    os.system("cd /glade/scratch/khalbert/ncarens/" + date + "/" + member + "; cp " + file + " ./;" + "echo n | gunzip /glade/scratch/khalbert/ncarens/" + date + "/" + member + "/" + file.split("/")[-1])
