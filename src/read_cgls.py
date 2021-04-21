#!/usr/bin/python
import os
import re
import datetime
import numpy as np 
import json
################################
# Fuctions to read CGLS data
################################
def get_cgls_locs(mapname="glb_15min"):
    # directory
    fname="/cluster/data6/menaka/CGLS/data/river/CGLS_alloc_"+mapname+".txt"
    #--
    nums=[]
    river=[]
    pname=[]
    lons =[]
    lats =[]
    xlist=[]
    ylist=[]
    egm08=[]
    egm96=[]
    #--
    f=open(fname,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        #print line
        num     = line[0]
        station = line[1]
        riv     = re.split("_",station)[1]
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[4])-1
        iy      = int(line[5])-1
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        #-----------------------
        nums.append(num)
        river.append(riv)
        pname.append(station)
        lons.append(lon)
        lats.append(lat)
        xlist.append(ix)
        ylist.append(iy)
        egm08.append(EGM08)
        egm96.append(EGM96)
    return nums,river,pname,lons,lats,xlist,ylist,egm08,egm96
################################
def cgls_continous_WSE(station,syear=2002,smon=10,sday=1,eyear=2020,emon=12,eday=31,egm08=0.0,egm96=0.0):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    #---
    time=int((end-start).days + 1) # time in days
    data=np.ones([time],np.float32)*-9999.0 # WSE in [m] # -9999.0 for no observations
    # read CGLS
    fname="/cluster/data6/menaka/CGLS/data/river/"+station+".json"
    with open(fname) as f:
        alldata    = json.load(f)
        cgls_data  = alldata["data"]
    #----------------------------
    for line in range(len(cgls_data)):
        date    = cgls_data[line]["datetime"]
        date    = re.split(" ",date)[0]
        date    = re.split("/",date)
        yyyy    = int(date[0])
        mm      = int(date[1])
        dd      = int(date[2])
        now     = datetime.date(yyyy,mm,dd)
        wse     = cgls_data[line]["water_surface_height_above_reference_datum"]
        if now > start and now < end:
            nowtime=int((now-start).days)
            data[nowtime]=wse+egm08-egm96
    return data
#####################################
def metadata():
    # directory 
    hydroweb="/cluster/data6/menaka/CGLS/data/river/CGLS_Station_list.txt"
    #--
    f=open(hydroweb,"r")
    lines=f.readlines()
    sta={}
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        id      = int(line[0])
        station = line[1]
        river   = line[2]
        basin   = line[3]
        country = line[4]
        sat     = line[10]
        startdt = line[11]
        enddt   = line[13]
        status  = line[15].split("\n")[0]
        #----
        sta[station]=[river,basin,country,sat,startdt,enddt,status]
    #----
    return sta
################################