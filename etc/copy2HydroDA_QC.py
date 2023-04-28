#!/opt/local/bin/python
# -*- coding: utf-8 -*-
"""
copying quality controlled allocation files to HydroDA
"""

import numpy as np
import sys
import os
import re
from multiprocessing import Pool
# from read_patchMS import upstream
############################################################
def slope(ix,iy,nextxy,uparea,elevtn,nxtdst,rivseq):
    nextX=nextxy[0]
    nextY=nextxy[1]
    slp1=0.0
    slp2=0.0
    if rivseq[iy,ix]>1 and nextX[iy,ix]>0:
        uXX, uYY = upstream(ix+1,iy+1,nextX.T,nextY.T,uparea.T)
        uXX = uXX - 1
        uYY = uYY - 1
        slp1=(elevtn[uYY,uXX]-elevtn[iy,ix])
        dXX = nextX[iy,ix] - 1
        dYY = nextY[iy,ix] - 1
        slp2=(elevtn[iy,ix]-elevtn[dYY,dXX])
    elif rivseq[iy,ix]==1:
        slp1=0.0
        dXX = nextX[iy,ix] - 1
        dYY = nextY[iy,ix] - 1
        slp2=(elevtn[iy,ix]-elevtn[dYY,dXX])
    elif nextX[iy,ix]<=0:
        uXX, uYY = upstream(ix+1,iy+1,nextX.T,nextY.T,uparea.T)
        uXX = uXX - 1
        uYY = uYY - 1
        slp1=(elevtn[uYY,uXX]-elevtn[iy,ix])
        slp2=0.0
    return slp1,slp2
############################################################
# dataname=sys.argv[1]
# mapname=sys.argv[2] #"glb_06min"
# glbname=sys.argv[3] #"glb_06min"
# CaMa_dir=sys.argv[4] #"/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# obstxt=sys.argv[5] #"./out/altimetry_"+mapname+"_test.txt"
# outtxt=sys.argv[6]
mapname="conus_06min"
# TAG="HydroWeb"
TAG="CGLS"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# intxt="../out/altimetry_"+mapname+"_20230327.txt"
intxt="../out/biased_removed_altimetry_"+mapname+"_20230406.txt"
outtxt="/cluster/data6/menaka/HydroDA/dat/"+TAG+"_alloc_"+mapname+"_DIR.txt"
############################################################
area_thr = 1.0e-20 #m2
slpe_thr = 1.0e20 #m
elev_thr = 1.0e20 #m
dist_thr = 1.0    #km
rmse_thr = 1.0    #m
bias_thr = 1.0    #m
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
nextxy = CaMa_dir+"/map/"+mapname+"/nextxy.bin"
rivwth = CaMa_dir+"/map/"+mapname+"/rivwth.bin"
rivhgt = CaMa_dir+"/map/"+mapname+"/rivhgt.bin"
rivlen = CaMa_dir+"/map/"+mapname+"/rivlen.bin"
elevtn = CaMa_dir+"/map/"+mapname+"/elevtn.bin"
uparea = CaMa_dir+"/map/"+mapname+"/uparea.bin"
lonlat = CaMa_dir+"/map/"+mapname+"/lonlat.bin"
nxtdst = CaMa_dir+"/map/"+mapname+"/nxtdst.bin"
rivseq = CaMa_dir+"/map/"+mapname+"/rivseq.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
rivseq = np.fromfile(rivseq,np.int32).reshape(ny,nx)
############################################################
#===============================================
# higher RMSE locations
fname="/cluster/data6/menaka/AltiMaP/out/rmse_VIC_BC_glb_06min_20211122.txt"
with open(fname,"r") as f_rmse:
	rmse=f_rmse.readlines()
#----
rmse_stations=[]
bias_stations=[]
for item in rmse[1::]:
    item    = list(filter(None, re.split(" ",item)))
    station = item[0]
    rmse0   = float(item[1])
    bias0   = float(item[2])
    if rmse0 > rmse_thr:
        rmse_stations.append(station)
        # print (rmse0, station)
    if bias0 > bias_thr:
        bias_stations.append(station)
        # print (bias0, station)
############################################################
#===============================================
# unreal observations
fname="../out/altimetry_"+mapname+"_20230327.txt"
with open(fname,"r") as f_unreal:
	unreal=f_unreal.readlines()
#----
unreal_stations=[]
for item in unreal:
    item    = list(filter(None, re.split(" ",item)))
    station = item[0]
    flag    = item[6]
    if flag >= 900:
        unreal_stations.append(station)
############################################################
#===========================================================
# open altimetry allocation file
with open(intxt,"r") as f:
	lines=f.readlines()
#=====================================
with open(outtxt,"w") as fout:
    fout.write("%13s%64s%12s%12s%8s%8s%12s%12s%12s%12s%17s\n"%("ID", "station", "lon", "lat", "ix", "iy", "elevation", "ele_diff", "EGM08", "EGM96", "satellite"))
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        # print (line)
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

        ################
        # condition for elevation difference
        ################
        eled=elevtn[iy-1,ix-1]-elev
        if abs(eled) > elev_thr:
            print ("elevation difference is too large: (>"+"%6.2f"%(elev_thr)+"m)", eled, "m")
            continue

        ################
        # condition for upstream catchment area
        ################
        if uparea[iy-1,ix-1] < area_thr:
            # print ("smaller river: (uparea <"+"%12.2f"%(area_thr*1e-6)+"$km^2$): ",uparea[iy-1,ix-1],station )
            continue
        
        ################
        # condition for slope
        ################
        # slp1, slp2 =slope(ix-1,iy-1,nextxy,uparea,elevtn,nxtdst,rivseq)
        # if slp1 > slpe_thr:
        #     print ("high slope:",station, slp1)
        #     continue

        ################
        # condition for distance to mouth
        ################
        if dist > dist_thr:
            print ("large distance to mouth: ",station, dist)
            continue

        ################
        # condition for unreal observations
        ################
        if station in unreal_stations:
            print ("unreal observations: ", station)
            continue

        ################
        # condition for higher rmse observations
        ################
        if station in rmse_stations:
            print ("higher RMSE virtual station ( >"+"%5.2f"%(rmse_thr)+"m): ", station)
            continue

        ################
        # condition for higher bias observations
        ################
        if station in bias_stations:
            print ("higher bias virtual station ( >"+"%5.2f"%(bias_thr)+"m): ", station)
            continue

        #========================
        linew="%13s%64s%12.2f%12.2f%8d%8d%12.2f%12.2f%12.2f%12.2f%17s\n"%(num,station,lon,lat,ix,iy,elev,eled,EGM08,EGM96,sat)
        # print (linew)
        # print (station, lon, lat, west, east, south, north)
        fout.write(linew)