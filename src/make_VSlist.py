#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Calculate the EGM08 and EGM96 values for each VS and write to the VS list
input: {datadir}/{dataname}_VS
identifier,river,basin,country,satellite,track_nb,start_date,end_date,latitude,longitude,status,validation
output: ./inp/{dataname}_Station_list.txt
ID | Station | River | Basin | Country | lon | lat | elevation | EGM08 | EGM96 | Satellite | Start Date | End Date | Status
"""
import os
import re
import sys
##########################
# get input
dataname=sys.argv[1]
datafile=sys.argv[2]
outdir=sys.argv[3]
##########################

##################
# read VS list
fname=datafile #+"/HydroWeb_VS"
with open(fname,"r") as fr:
    lines=fr.readlines()
#-
ID_list=[]
station_list=[]
river_list=[]
basin_list=[]
country_list=[]
sat_list=[]
start_list=[]
end_list=[]
status_list=[]
ele_list=[]
lat_list=[]
lon_list=[]
#============================
#  write VS lon lat INPUT.DAT 
#============================
with open("./INPUT.DAT","w") as fw:
    for line in lines[1::]:
        #print line
        line   = re.split(",",line)
        #line   = filter(None,line)
        station= line[0].strip()
        river  = line[1].strip()
        basin  = line[2].strip()
        country= line[3].strip()
        if country == "":
            country="Country"
        sat    = line[4].strip()
        lat    = float(line[8])
        lon    = float(line[9])
        startdt= line[6].split(' ')[0]
        enddt  = line[7].split(' ')[0]
        status = line[10].split('\n')[0]
        # lon lat
        linew = "%7.4f   %7.4f\n"%(lat,lon)
        print (linew)
        fw.write(linew)
        #===========================
        #--get ID & elevation
        #===========================
        if dataname == "HydroWeb":
            iname=datadir+"/hydroprd_"+station+".txt"
            with open(iname,"r") as f_hyd:
                l_hyd=f_hyd.readlines()
            ID=int(l_hyd[2].split("::")[-1])
            ele=float(l_hyd[16].split("::")[-1])
        else:
            print (station)
            ID=int(station[-4::])
            ele=-99.0
        #===========================
        # append
        #===========================
        station_list.append(station)
        ID_list.append(ID)
        river_list.append(river)
        basin_list.append(basin)
        country_list.append(country)
        sat_list.append(sat)
        start_list.append(startdt)
        end_list.append(enddt)
        status_list.append(status)
        ele_list.append(ele)
        lat_list.append(lat)
        lon_list.append(lon)
#+++++++++++++++++++++++
########################
#   get EGM08 values
########################
# run intpt_EGM08 
os.system("./src/intpt_EGM08")
#========================
# read output file
#========================
fr=open("./OUTPUT_EGM08.DAT","r")
liner=fr.readlines()
fr.close()
EGM08_list=[]
for i,line in enumerate(liner):
    line  =  re.split(" ",line)
    line  = filter(None,line)
    EGM08 = float(line[2])
    EGM08_list.append(EGM08)
    print (EGM08)
########################
#   get EGM96 values
########################
# run intpt_EGM96 
os.system("./src/intpt_EGM96")
#========================
# read output.dat
#========================
fr=open("./OUTPUT_EGM96.DAT","r")
liner=fr.readlines()
fr.close()
EGM96_list=[]
for i,line in enumerate(liner):
    line  =  re.split(" ",line)
    line  = filter(None,line)
    EGM96 = float(line[2])
    EGM96_list.append(EGM96)
    print (EGM96)
#==========================================
# write VS list
#==========================================
fname=outdir+"/"+dataname+"Station_list.txt"
with open(fname,"w") as fww:
    header="%15s%62s%32s%32s%32s%10s%10s%10s%10s%10s%17s%17s%17s%13s\n"%("ID","Station","River","Basin","Country","lon","lat","elevation","EGM08","EGM96","Satellite","Start Date","End Date","Status")
    fww.write(header)
    pnum=len(ID_list)
    print (pnum)
    for i in range(pnum):
        #print i, ID_list[i],station_list[i],river_list[i],basin_list[i],country_list[i],lon_list[i],lat_list[i],ele_list[i],sat_list[i],start_list[i],end_list[i],status_list[i]
        country='-'.join(country_list[i].split())
        linew="%015d%62s%32s%32s%32s%10.4f%10.4f%10.4f%10.4f%10.4f%17s%17s%17s%13s\n"%(ID_list[i],station_list[i],river_list[i],basin_list[i],country,lon_list[i],lat_list[i],ele_list[i],EGM08_list[i],EGM96_list[i],sat_list[i],start_list[i],end_list[i],status_list[i])
        print (linew)
        fww.write(linew)
#==========================================
os.system("rm ./INPUT.DAT")
os.system("rm ./OUTPUT_EGM08.DAT")
os.system("rm ./OUTPUT_EGM96.DAT")
print ("---> "+dataname+"Station_list.txt"+" done.")