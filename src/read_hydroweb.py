#!/usr/bin/python
import os
import re
import datetime
import numpy as np 
################################
# Fuctions to read HydroWeb data
################################
def hydroweb_river_name(mapname="glb_15min"):
    # directory
    hydroweb="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    #--
    river=[]
    #--
    f=open(hydroweb,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        print line
        station = line[1]
        riv     = re.split("_",station)[1]
        river.append(riv)
    return river
################################
def get_hydroweb(mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    #--
    river=[]
    pname=[]
    xlist=[]
    ylist=[]
    egm_d=[]
    #--
    f=open(fname,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        print line
        station = line[1]
        riv     = re.split("_",station)[1]
        ix      = int(line[4])-1
        iy      = int(line[5])-1
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        river.append(riv)
        pname.append(station)
        xlist.append(ix)
        ylist.append(iy)
        egm_d.append(EGM96-EGM08)
    return river,pname,xlist,ylist,egm_d
##################################
def get_hydroweb_loc(rivername,mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    # directory
    # fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    #--
    #river=[]
    pname=[]
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
        station = line[1]
        riv     = re.split("_",station)[1]
        ix      = int(line[4])-1
        iy      = int(line[5])-1
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        if rivername==riv:
            #river.append(riv)
            pname.append(station)
            xlist.append(ix)
            ylist.append(iy)
            egm08.append(EGM08)
            egm96.append(EGM96)
    return pname,xlist,ylist,egm08,egm96
##################################
def get_hydroweb_locs(mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    # directory
    # fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    # fname="./out/altimetry_"+mapname+"_test.txt"
    #--
    nums=[]
    river=[]
    pname=[]
    lons =[]
    lats =[]
    xlist=[]
    ylist=[]
    leled=[]
    egm08=[]
    egm96=[]
    llsat=[]
    ldtom=[]
    lflag=[]
    #--
    f=open(fname,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        #print line
        num     = line[0]
        station = line[1]
        riv     = re.split("_",station)[1]
        lon     = float(line[3])
        lat     = float(line[4])
        ix      = int(line[5])-1
        iy      = int(line[6])-1
        eled    = float(line[7])
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        sat     = line[10].strip()
        dist    = float(line[11])
        flag    = int(line[12])
        #-----------------------
        nums.append(num)
        river.append(riv)
        pname.append(station)
        lons.append(lon)
        lats.append(lat)
        xlist.append(ix)
        ylist.append(iy)
        leled.append(eled)
        egm08.append(EGM08)
        egm96.append(EGM96)
        llsat.append(sat)
        ldtom.append(dist)
        lflag.append(flag)
    return nums,river,pname,lons,lats,xlist,ylist,egm08,egm96,llsat,leled,ldtom,lflag
################################
def HydroWeb_WSE(station,syear,eyear,smon=1,emon=12,sday=1,eday=31,egm08=0.0,egm96=0.0,fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    # read hydroweb
    #station="R_con_con_env_0429_01"
    #satellite=station.split("_")[2]
    #fname="/home/yamadai/data/Altimetry/HydroWeb_LEGOS/River/R_"+station
    # fname="/cluster/data6/menaka/HydroWeb/data/hydroprd_"+station+".txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    head=33
    #--
    time=[] # time in days
    data=[] # WSE in [m]
    for line in lines[head::]:
        if line[0][0] == "#":
            continue
        line = re.split(" ",line)
        date = line[0]
        date = re.split("-",date)
        yyyy = int(date[0])
        mm   = int(date[1])
        dd   = int(date[2])
        if float(line[2]) >= 9999.0:
            continue
        wse  = float(line[2])+egm08-egm96
        #print yyyy, mm, dd, wse
        now  = datetime.date(yyyy,mm,dd)
        if now < start and now > end:
            continue
        data.append(wse)
        lag  = int((now-start).days)
        time.append(lag)
    return time, data
#####################################
def HydroWeb_continous_WSE(station,syear=2002,smon=10,sday=1,eyear=2020,emon=12,eday=31,egm08=0.0,egm96=0.0,dir0="/cluster/data6/menaka/HydroWeb/data"):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    # read hydroweb
    #station="R_con_con_env_0429_01"
    #satellite=station.split("_")[2]
    #fname="/home/yamadai/data/Altimetry/HydroWeb_LEGOS/River/R_"+station
    # fname="/cluster/data6/menaka/HydroWeb/data/"hydroprd_"+station+".txt"
    fname=dir0+"/hydroprd_"+station+".txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    head=33
    #--
    time=int((end-start).days + 1) # time in days
    data=np.ones([time],np.float32)*-9999.0 # WSE in [m] # -9999.0 for no observations
    for line in lines[head::]:
        if line[0][0] == "#":
            continue
        line = re.split(" ",line)
        date = line[0]
        date = re.split("-",date)
        yyyy = int(date[0])
        mm   = int(date[1])
        dd   = int(date[2])
        if float(line[2]) >= 9999.0:
            continue
        wse  = float(line[2])
        #print yyyy, mm, dd, wse
        now  = datetime.date(yyyy,mm,dd)
        #print yyyy, mm, dd, int((now-start).days), wse
        if now > start and now < end:
            nowtime=int((now-start).days)
            data[nowtime]=wse+egm08-egm96
    return data
#####################################
def altimetry(name,mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    # directory
    hydroweb="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    #--
    f=open(fname,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        station = line[1]
        if station==name:
            alti= float(line[6])
    return alti
################################
def metadata():
    # directory 
    hydroweb="/cluster/data6/menaka/HydroWeb/HydroWebStation_list.txt"
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
        try:
            status  = line[15].split("\n")[0]
        except:
            status = "N/A"
        #----
        sta[station]=[river,basin,country,sat,startdt,enddt,status]
    #----
    return sta
################################