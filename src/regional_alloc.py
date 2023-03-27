#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import sys
import re
import numpy as np
'''
Covert allocated locations for regional map [e.g. amz_06min etc.]
'''
#===================
dataname=sys.argv[1]
mapname=sys.argv[2] #"glb_06min"
glbname=sys.argv[3] #"glb_06min"
CaMa_dir=sys.argv[4] #"/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
obstxt=sys.argv[5] #"./out/altimetry_"+mapname+"_test.txt"
outtxt=sys.argv[6] #
# odir=sys.argv[2] #"/cluster/data6/menaka/Altimetry/results"
# restag=sys.argv[5] #"3sec"
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
# global map
fname=CaMa_dir+"/map/"+glbname+"/params.txt"
with open(fname,"r") as fmap:
    lines=fmap.readlines()
#-------
# nx     = int(filter(None, re.split(" ",lines[0]))[0])
# ny     = int(filter(None, re.split(" ",lines[1]))[0])
# gsize  = float(filter(None, re.split(" ",lines[3]))[0])
lon_ori  = float(filter(None, re.split(" ",lines[4]))[0])
lon_end  = float(filter(None, re.split(" ",lines[5]))[0])
lat_end  = float(filter(None, re.split(" ",lines[6]))[0])
lat_ori  = float(filter(None, re.split(" ",lines[7]))[0])
############################################################
dx=int( (west-lon_ori) /gsize +0.001 ) ##  add 0.001 to avoid rounding error
dy=int( (lat_ori-north)/gsize +0.001 )
print ("dx: ",dx,"dy: ",dy)
############################################################
#-------------------------------------------
# fname="./out/altimetry_"+mapname+"_test.txt"
# fname="./out/altimetry_"+mapname+"_20210518.txt"
fname=obstxt
#--
with open(fname,"r") as f:
    lines=f.readlines()
#=====================================
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
    ix      = int(line[16])-dx
    iy      = int(line[17])-dy
    EGM08   = float(line[18])
    EGM96   = float(line[19])

    # num     = line[0]
    # station = line[1]
    # line2   = re.split("_",station)
    # riv     = line2[1]
    # stream  = line2[2]
    # data    = line[2]
    # lon     = float(line[3])
    # lat     = float(line[4])
    # ix      = int(line[5])-dx
    # iy      = int(line[6])-dy
    # eled    = float(line[7])
    # EGM08   = float(line[8])
    # EGM96   = float(line[9])
    # sat     = line[10].strip()
    # dist    = float(line[11])
    # flag    = int(line[12])
    # kx1     = int(line[13])
    # ky1     = int(line[14])
    # kx2     = int(line[15])
    # ky2     = int(line[16])
    # dist1   = float(line[17])
    # dist2   = float(line[18])
    #======================
    if ix < 1:
        continue
    if ix > nx:
        continue 
    if iy < 1:
        continue
    if iy > ny:
        continue 
    if lon >= west and lon <= east and lat >= south and lat <= north:
        print ("%13s%64s%12s%12.2f%12.2f%17s%6d%12.2f%15.2f%10d%8d%8d%8d%14.2f%12.2f%12.2f%10d%8d%12.2f%10.2f")%(num,station,dataname,lon,lat,sat,flag,elev,dist,kx1,ky1,kx2,ky2,dist1,dist2,rivwth,ix,iy,EGM08,EGM96)
        # linew="%30s%67s%12s%12.2f%12.2f%8d%8d%12.2f%12.2f%12.2f%17s%12.2f%6d%10d%10d%10d%10d%12.2f%12.2f\n"%(num,station,data,lon,lat,ix,iy,eled,EGM08,EGM96,sat,dist,flag,kx1,ky1,kx2,ky2,dist1,dist2)
        # print (linew)
# # outtxt="./out/altimetry_"+mapname+"_test.txt"
# fname=outtxt
# #--
# with open(fname,"w") as fout:
#     header = "%30s%67s%12s%12s%10s%10s%8s%12s%12s%12s%17s%15s%6s%10s%8s%8s%8s%12s%12s\n"%("ID","station","dataname","lon","lat","ix","iy","elevation","EGM08","EGM96","satellite","dist_to_mouth","flag","kx1","ky1","kx2","ky2","dist1","dist2")
#     fout.write(header)
#     for line in lines[1::]:
#         line    = filter(None,re.split(" ",line))
#         # print line
#         num     = line[0]
#         station = line[1]
#         line2   = re.split("_",station)
#         riv     = line2[1]
#         stream  = line2[2]
#         data    = line[2]
#         lon     = float(line[3])
#         lat     = float(line[4])
#         ix      = int(line[5])-dx
#         iy      = int(line[6])-dy
#         eled    = float(line[7])
#         EGM08   = float(line[8])
#         EGM96   = float(line[9])
#         sat     = line[10].strip()
#         dist    = float(line[11])
#         flag    = int(line[12])
#         kx1     = int(line[13])
#         ky1     = int(line[14])
#         kx2     = int(line[15])
#         ky2     = int(line[16])
#         dist1   = float(line[17])
#         dist2   = float(line[18])
#         #======================
#         if ix < 1:
#             continue
#         if ix > nx:
#             continue 
#         if iy < 1:
#             continue
#         if iy > ny:
#             continue 
#         if lon >= west and lon <= east and lat >= south and lat <= north:
#             linew="%30s%67s%12s%12.2f%12.2f%8d%8d%12.2f%12.2f%12.2f%17s%12.2f%6d%10d%10d%10d%10d%12.2f%12.2f\n"%(num,station,data,lon,lat,ix,iy,eled,EGM08,EGM96,sat,dist,flag,kx1,ky1,kx2,ky2,dist1,dist2)
#             # print (linew)
#             print (station, lon, lat, west, east, south, north)
#             fout.write(linew)