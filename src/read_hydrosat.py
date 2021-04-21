#!/usr/bin/python
import os
import re
import datetime
import numpy as np 
################################
# Fuctions to read HydroSat data
################################
def hydrosat_river_name(mapname="glb_15min"):
    # directory
    hydrosat="/cluster/data6/menaka/HydroSat/HydroSat_alloc_"+mapname+".txt"
    #--
    river=[]
    #--
    f=open(hydrosat,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        print line
        station = line[1]
        riv     = re.split("_",station)[1]
        river.append(riv)
    return river
################################
def get_hydrosat(mapname="glb_15min"):
    # directory
    hydrosat="/cluster/data6/menaka/HydroSat/HydroSat_alloc_"+mapname+".txt"
    #--
    river=[]
    pname=[]
    xlist=[]
    ylist=[]
    egm_d=[]
    #--
    f=open(hydrosat,"r")
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
def get_hydroweb_loc(rivername,mapname="glb_15min"):
    # directory
    hydrosat="/cluster/data6/menaka/HydroSat/HydroSat_alloc_"+mapname+".txt"
    #--
    #river=[]
    pname=[]
    xlist=[]
    ylist=[]
    egm08=[]
    egm96=[]
    #--
    f=open(hydrosat,"r")
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
def get_hydrosat_locs(mapname="glb_15min"):
    # directory
    fname="/cluster/data6/menaka/HydroSat/HydroSat_alloc_"+mapname+".txt"
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
        num     = line[0].strip()
        station = line[1].strip()
        riv     = line[2].strip()
        lon     = float(line[2])
        lat     = float(line[3])
        ix      = int(line[4])-1
        iy      = int(line[5])-1
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        #print num, station, riv
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
def HydroSat_WSE(station,syear,eyear,smon=1,emon=12,sday=1,eday=31,egm08=0.0,egm96=0.0):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    # read HydroSat
    id=get_id(station)
    fname="/cluster/data6/menaka/HydroSat/data/%d.txt"%(id)
    with open(fname) as f:
        lines= [line.strip() for line in f if line.strip()]
    # f=open(fname,"r")
    # lines=f.readlines()
    # f.close()
    head=35
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
def HydroSat_continous_WSE(station,syear=2002,smon=10,sday=1,eyear=2020,emon=12,eday=31,egm08=0.0,egm96=0.0):
    #==================================
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    # read HydroSat
    id=get_id(station)
    #id=station
    fname="/cluster/data6/menaka/HydroSat/data/%d.txt"%(id)
    with open(fname) as f:
        lines= [line.strip() for line in f if line.strip()]
    # f=open(fname,"r")
    # lines=f.readlines()
    # f.close()
    head=35
    #--
    time=int((end-start).days + 1) # time in days
    data=np.ones([time],np.float32)*-9999.0 # WSE in [m] # -9999.0 for no observations
    for line in lines[head::]:
        if line[0][0] == "#":
            continue
        #print line
        line = re.split(",",line)
        yyyy = int(line[0])
        mm   = int(line[1])
        dd   = int(line[2])
        if line[3].strip()=="NaN":
            continue
        wse  = float(line[3])
        #print yyyy, mm, dd, wse
        now  = datetime.date(yyyy,mm,dd)
        #print yyyy, mm, dd, int((now-start).days), wse
        if now > start and now < end:
            nowtime=int((now-start).days)
            data[nowtime]=wse+egm08-egm96
    return data
#####################################
def altimetry(name,mapname="glb_15min"):
    # directory
    hydrosat="/cluster/data6/menaka/HydroSat/HydroSat_alloc_"+mapname+".txt"
    #--
    f=open(hydrosat,"r")
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
    hydrosat="/cluster/data6/menaka/HydroSat/HydroSatStation_list.txt"
    #--
    f=open(hydrosat,"r")
    lines= [line.strip() for line in f if line.strip()]
    #lines=f.readlines()
    sta={}
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        id      = int(line[0])
        station = line[1].strip()
        river   = line[2].strip()
        basin   = line[3].strip()
        #print id,station,river
        country = line[4]
        sat     = line[10]
        startdt = line[11]
        enddt   = line[12]
        status  = line[13].split("\n")[0]
        #----
        sta[station]=[river,basin,country,sat,startdt,enddt,status]
    #----
    return sta
################################
def get_id(name):
    # directory 
    hydrosat="/cluster/data6/menaka/HydroSat/HydroSatStation_list.txt"
    #--
    f=open(hydrosat,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        id      = int(line[0])
        station = line[1].strip()
        #print id, station, name
        if name.strip()==station.strip():
            #print "*******", id, station, name
            HydroSat_id = id
            break
    return HydroSat_id
################################