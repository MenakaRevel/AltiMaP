#!/usr/bin/python
import os
import re
import datetime
import numpy as np 
import json
################################
# Fuctions to read GRRATS data
################################
def GRRATS_WSE(station,syear=2002,eyear=2020,smon=1,emon=12,sday=1,eday=31,egm08=0.0,egm96=0.0):
    #----
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)
    # read GRRATS
    fname="/cluster/data6/menaka/GRRATS/txtdata/"+station+".txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    head=19
    #--
    time=[] # time in days
    data=[] # WSE in [m]
    for line in lines[head::]:
        if line[0][0] == "#":
            continue
        line = filter(None,re.split(" ",line))
        # print line
        date = line[0]
        date = re.split("-",date)
        yyyy = int(date[0])
        mm   = int(date[1])
        dd   = int(date[2])
        if line[3].strip()=="nan":
            continue
        wse  = float(line[3])+egm08-egm96
        #print yyyy, mm, dd, wse
        now  = datetime.date(yyyy,mm,dd)
        if now < start and now > end:
            continue
        data.append(wse)
        lag  = int((now-start).days)
        time.append(lag)
    return time, data
#####################################