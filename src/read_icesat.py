#!/usr/bin/python
import os
import re
import datetime
import numpy as np 
################################
# Fuctions to read ICESat data
################################
def get_icesat_locs(mapname="glb_15min"):
    # directory
    fname="/cluster/data6/menaka/ICEsat/ICESat_"+mapname+"/ICESat_list.txt"
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
        num     = 10000*int(line[0])+int(line[1])
        station = "%04d%04d"%(int(line[0]),int(line[1]))
        riv     = "River" #re.split("_",station)[1]
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[0])-1
        iy      = int(line[1])-1
        EGM08   = 0.0 #float(line[8])
        EGM96   = 0.0 #float(line[9])
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
#####################################
def ICESat_continous_WSE(station,syear=2003,smon=1,sday=1,eyear=2020,emon=12,eday=31,egm08=0.0,egm96=0.0):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    # read hydroweb
    #station="R_con_con_env_0429_01"
    #satellite=station.split("_")[2]
    #fname="/home/yamadai/data/Altimetry/HydroWeb_LEGOS/River/R_"+station
    fname="/cluster/data6/menaka/ICEsat/ICESat_glb_06min/"+station+".txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    head=12
    #--
    time=int((end-start).days + 1) # time in days
    data=np.ones([time],np.float32)*-9999.0 # WSE in [m] # -9999.0 for no observations
    for line in lines[head::]:
        if line[0][0] == "#":
            continue
        line = filter(None,re.split(" ",line))
        #print line
        date = line[0]
        date = re.split("-",date)
        yyyy = int(date[0])
        mm   = int(date[1])
        dd   = int(date[2])
        wse  = float(line[2])
        #print yyyy, mm, dd, wse
        now  = datetime.date(yyyy,mm,dd)
        #print yyyy, mm, dd, int((now-start).days), wse
        if now > start and now < end:
            nowtime=int((now-start).days)
            data[nowtime]=wse+egm08-egm96
    return data
#####################################
def metadata():
    # directory 
    hydroweb="/cluster/data6/menaka/ICEsat/ICESat_glb_06min/ICESat_list.txt"
    #--
    f=open(hydroweb,"r")
    lines=f.readlines()
    sta={}
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        id      = 10000*int(line[0])+int(line[1]) #int(line[0])
        station = "%04d%04d"%(int(line[0]),int(line[1])) #line[1]
        river   = "River" #line[2]
        basin   = "River" #line[3]
        country = "River" #line[4]
        sat     = "ICESat 1" #line[10]
        startdt = "2003-1-1"#line[11]
        enddt   = "2009-12-31"#line[13]
        status  = "ICESat 1"#line[15].split("\n")[0]
        #----
        sta[station]=[river,basin,country,sat,startdt,enddt,status]
    #----
    return sta
################################