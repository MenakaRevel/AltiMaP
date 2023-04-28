#!/opt/local/bin/python
# -*- coding: utf-8 -*-
"""
copying original allocation files to HydroDA
"""

import numpy as np
import sys
import os
import re
from multiprocessing import Pool
############################################################
# dataname=sys.argv[1]
# mapname=sys.argv[2] #"glb_06min"
# glbname=sys.argv[3] #"glb_06min"
# CaMa_dir=sys.argv[4] #"/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# obstxt=sys.argv[5] #"./out/altimetry_"+mapname+"_test.txt"
# outtxt=sys.argv[6]
mapname="conus_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# intxt="../out/altimetry_"+mapname+"_20230327.txt"
# intxt="../out/altimetry_"+mapname+"_20230404.txt"
intxt="../out/biased_removed_altimetry_"+mapname+"_20230406.txt"
# out/biased_removed_altimetry_glb_06min_20230405.txt
outtxt="/cluster/data6/menaka/HydroDA/dat/CGLS_alloc_"+mapname+"_org.txt"
############################################################
# regional map
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as fmap:
    lines=fmap.readlines()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])
############################################################
elevtn = CaMa_dir+"/map/"+mapname+"/elevtn.bin"
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
############################################################
with open(intxt,"r") as f:
    lines=f.readlines()
#=====================================
with open(outtxt,"w") as fout:
    fout.write("%13s%64s%12s%12s%8s%8s%12s%12s%12s%12s%17s\n"%("ID", "station", "lon", "lat", "ix", "iy", "elevation", "ele_diff", "EGM08", "EGM96", "satellite"))
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        # print line
        num     = line[0]
        station = line[1]
        line2   = re.split("_",station)
        riv     = line2[1]
        stream  = line2[2]
        dataname= line[2]
        lon     = float(line[3])
        lat     = float(line[4])
        sat     = line[5]
        #--
        flag    = int(line[6])
        elev    = float(line[7])
        dist    = float(line[8])
        kx1     = int(line[9])
        ky1     = int(line[10])
        kx2     = int(line[11])
        ky2     = int(line[12])
        dist1   = float(line[13])
        dist2   = float(line[14])
        rivwth  = float(line[15])
        #--
        ix      = int(line[16])
        iy      = int(line[17])
        EGM08   = float(line[18])
        EGM96   = float(line[19])
        #========================
        eled=elevtn[iy-1,ix-1]-elev
        linew="%13s%64s%12.2f%12.2f%8d%8d%12.2f%12.2f%12.2f%12.2f%17s\n"%(num,station,lon,lat,ix,iy,elev,eled,EGM08,EGM96,sat)
        # print (linew)
        print (station, lon, lat, west, east, south, north)
        fout.write(linew)