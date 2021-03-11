#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
# import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import matplotlib.patches as mpatches
import re
import os
import sys
import errno
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes

import read_hydroweb as hweb
import read_cgls as cgls
import read_hydrosat as hsat
import read_icesat as isat
import read_grrats as grt
#=============================
def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
#=============================
def get_data(station, tag):
    time=-9999.0
    data=-9999.0
    if tag=="HydroWeb":
        time,data=hweb.HydroWeb_WSE(station,syear=1992)
    elif tag=="HydroSat":
        time,data=hsat.HydroSat_WSE(station,syear=1992)
    elif tag=="ICESat":
        time,data=isat.ICESat_WSE(station,syear=1992)
    elif tag=="CGLS":
        time,data=cgls.CGLS_WSE(station,syear=1992)
    elif tag=="GRRATS":
        time,data=grt.GRRATS_WSE(station,syear=1992)
    return time, data
#=============================
def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5) / multiplier
#=============================
def round_half_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n*multiplier - 0.5) / multiplier
#=============================
def vec_par(inlist):
    LEVEL = inlist[0]
    ixiy  = inlist[1]
    ax=None
    sup=2
    w=0.02
    width=0.5
    ax=ax or plt.gca()
    txt="%s_%s.txt"%(ixiy,LEVEL)
    os.system("./src/print_rivvec "+ixiy+".txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    with open(txt,"r") as f:
        lines = f.readlines()
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        if lon1-lon2 > 180.0:
            print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            print (lon1,lon2)
            lon2=-180.0

        colorVal="xkcd:azure"
        # print lon1,lon2,lat1,lat2
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal)
#=============================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None,alpha=1):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha) #transform=ccrs.PlateCarree(),
#=============================
syear=1992
eyear=2020
start=datetime.date(syear,1,1)
end=datetime.date(eyear,12,31)
days=(end-start).days + 1
#=============================
idnames={}
pnames={}
dataname={}
listnum={}
xxlist={}
yylist={}
satlist={}
lonlist={}
latlist={}
#=============================
# mapname="glb_01min"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# print sys.argv
mapname=sys.argv[1]
CaMa_dir=sys.argv[2]
ncpus=int(sys.argv[3])
fname="./out/altimetry_"+mapname+".txt"
with open(fname, "r") as f:
    lines=f.readlines()
    for line in lines[1::]:
        line = filter(None,re.split(" ", line))
        # print line
        Id   = line[0].strip()
        name = line[1].strip()
        dname= line[2].strip()
        lon  = float(line[3])
        lat  = float(line[4])
        ix   = int(line[5]) 
        iy   = int(line[6])
        eled = float(line[7])
        egm08= float(line[8])
        egm96= float(line[9])
        sat  = line[10].strip()
        ixiy="%05d%05d"%(ix,iy)
        #print ixiy
        if ixiy in pnames.keys():
            pnames[ixiy].append(name)
            dataname[ixiy].append(dname)
            xxlist[ixiy].append(ix)
            yylist[ixiy].append(iy)
            lonlist[ixiy].append(lon)
            latlist[ixiy].append(lat)
            satlist[ixiy].append(sat.upper())
        else:
            pnames[ixiy]=[name]
            dataname[ixiy]=[dname]
            xxlist[ixiy]=[ix]
            yylist[ixiy]=[iy]
            lonlist[ixiy]=[lon]
            latlist[ixiy]=[lat]
            satlist[ixiy]=[sat.upper()]
#=============================
# Read the CMF variables
if mapname == 'glb_15min':
    nx      = 1440
    ny      = 720
elif mapname == 'glb_06min':
    nx      = 3600
    ny      = 1800
elif mapname == 'glb_01min':
    nx      = 21600
    ny      = 10800
#=============================
nextxy = CaMa_dir+"/map/"+mapname+"/nextxy.bin"
rivwth = CaMa_dir+"/map/"+mapname+"/rivwth.bin"
rivhgt = CaMa_dir+"/map/"+mapname+"/rivhgt.bin"
rivlen = CaMa_dir+"/map/"+mapname+"/rivlen.bin"
elevtn = CaMa_dir+"/map/"+mapname+"/elevtn.bin"
lonlat = CaMa_dir+"/map/"+mapname+"/lonlat.bin"
uparea = CaMa_dir+"/map/"+mapname+"/uparea.bin"
nxtdst = CaMa_dir+"/map/"+mapname+"/nxtdst.bin"
rivseq = CaMa_dir+"/map/"+mapname+"/rivseq.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
# rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
# rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
# rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
rivseq = np.fromfile(rivseq,np.int32).reshape(ny,nx)
#=====================================
markers={"HydroWeb":"o","CGLS":"s","ICESat":"^","HydroSat":"X","GRRATS":"D"}
colors={"HydroWeb":"xkcd:reddy brown","CGLS":"xkcd:dark pink","ICESat":"xkcd:pinkish","HydroSat":"xkcd:light urple","GRRATS":"xkcd:tangerine"} 
#=============================
maps = ['ESRI_Imagery_World_2D',    # 0
        'ESRI_StreetMap_World_2D',  # 1
        'NatGeo_World_Map',         # 2
        'NGS_Topo_US_2D',           # 3
        'Ocean_Basemap',            # 4
        'USA_Topo_Maps',            # 5
        'World_Imagery',            # 6
        'World_Physical_Map',       # 7
        'World_Shaded_Relief',      # 8
        'World_Street_Map',         # 9
        'World_Terrain_Base',       # 10
        'World_Topo_Map'            # 11
        ]
#=============================
# mkdir("./fig")
# mkdir("./fig/river_network")
#================================
#=== function for writing pdf ===
#================================
def write_pdf(inputlist):
    num   = inputlist[0]
    start = int(inputlist[1])
    last  = int(inputlist[2])
    pdfname="./fig/river_network/WSE_observation_river_network_"+num+".pdf"
    with PdfPages(pdfname) as pdf:
        for ixiy in pnames.keys()[start:last]:
            unilist=np.unique(np.array(dataname[ixiy]))
            pnum=np.shape(np.unique(np.array(dataname[ixiy])))[0]
            # if pnum<2:
            #     continue
            alpha=1
            sup=2
            w=0.02
            width=0.5
            hgt=11.69
            wdt=8.27
            fig=plt.figure(figsize=(wdt, hgt))
            #plt.title(pname[point][0],fontsize=12)
            G    = gridspec.GridSpec(3,2)
            ax00 = fig.add_subplot(G[0:2,0::])
            #ax01 = fig.add_subplot(G[0,1])
            ax11 = fig.add_subplot(G[2,:])
            # get the dimesion of the map
            dec=2
            val=0.10
            lllat = round_half_down(min(latlist[ixiy])-val,dec)
            urlat = round_half_up(max(latlist[ixiy])+val,dec)
            lllon = round_half_down(min(lonlist[ixiy])-val,dec)
            urlon = round_half_up(max(lonlist[ixiy])+val,dec)
            if abs(lllat-urlat) < val:
                urlat=round_half_up(urlat+val,dec)
                lllat=round_half_down(lllat-val,dec)
            if abs(lllon-urlon) < val:
                urlon=round_half_up(urlon+val,dec)
                lllon=round_half_down(lllon-val,dec)
            print (lllat, lllon, urlat, urlon)
            M = Basemap(resolution='h', projection='cyl',llcrnrlon=lllon, llcrnrlat=lllat, \
                urcrnrlon=urlon, urcrnrlat=urlat, ax=ax00)
            try:
                M.arcgisimage(service=maps[0], xpixels=1500, verbose=False)
            except:
                # Draw some map elements on the map
                M.drawcoastlines()
                M.drawstates()
                M.drawcountries()
                M.drawrivers(color='blue')
            M.drawparallels([lllat,urlat], labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
            M.drawmeridians([lllon,urlon], labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
            #####
            box="%f %f %f %f"%(lllon,urlon,urlat,lllat) 
            # print box
            os.system("./src/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > ./tmp/"+ixiy+".txt") 
            inlist=[]
            for LEVEL in np.arange(1,10+1,1):
                txt="./tmp/%s_%02d.txt"%(ixiy,LEVEL)
                os.system("./src/print_rivvec ./tmp/"+ixiy+".txt 1 "+str(LEVEL)+" > "+txt)
                width=(float(LEVEL)**sup)*w
                with open(txt,"r") as f:
                    lines = f.readlines()
                #---
                for line in lines:
                    line    = filter(None, re.split(" ",line))
                    lon1 = float(line[0])
                    lat1 = float(line[1])
                    lon2 = float(line[3])
                    lat2 = float(line[4])

                    if lon1-lon2 > 180.0:
                        print (lon1, lon2)
                        lon2=180.0
                    elif lon2-lon1> 180.0:
                        print (lon1,lon2)
                        lon2=-180.0

                    colorVal="xkcd:azure"
                    # print lon1,lon2,lat1,lat2
                    ax00.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
            # map(vec_par,inlist)
            os.system("rm -r ./tmp/"+ixiy+"*.txt")
            for i in np.arange(pnum):
                TAG=unilist[i] 
                repeatlist=np.where(np.array(dataname[ixiy])==TAG)[0]
                for j in repeatlist:
                    pname=pnames[ixiy][j]
                    lon0=lonlist[ixiy][j]
                    lat0=latlist[ixiy][j]
                    ix=int(ixiy[0:5])
                    iy=int(ixiy[5::])
                    lon=lonlat[0,iy-1,ix-1]
                    lat=lonlat[1,iy-1,ix-1]
                    # M.scatter(lon,lat,c=colors[TAG],s=20,marker=markers[TAG],zorder=110)
                    ax11.plot(lon0,lat0,color=colors[TAG],marker=markers[TAG],markersize=5)
                    ax11.plot(lon ,lat ,color=colors[TAG],marker=markers[TAG],fillstyle="none",markersize=5)
                    ax11.text(-0.1,1.1-0.1*j,pnames[ixiy][j],va="center",ha="left",transform=ax11.transAxes,fontsize=10)
            ax11.set_axis_off()
            ax11.spines['top'].set_visible(False)
            ax11.spines['right'].set_visible(False)
            ax11.spines['left'].set_visible(False)
            ax11.spines['bottom'].set_visible(False)
            ######################################
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
            print ("============================")
        # set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Title'] = 'Comparison of HydroWeb, CGLS, HydroSat, ICESat along the river'
        d['Author'] = 'Menaka Revel'
        d['Subject'] = 'Comparison of altimetry observations'
        d['Keywords'] = 'HydroWeb, CGLS, HydroSat, ICESat, GRRATS'
        d['CreationDate'] = datetime.datetime(2021, 1, 25)
        d['ModDate'] = datetime.datetime.today()
    return 0
#==============================================
pages=100.0
nums = int(len(pnames.keys())/pages) + 1
lastitem = len(pnames.keys())
inputlist=[]
for num in np.arange(0,nums):
    start = num*pages
    last  = (num + 1)*pages
    if last >= len(pnames.keys()):
        last = lastitem
    numch   = "%03d"%(num+1)
    startch = "%d"%(start)
    lastch  = "%d"%(last)
    inputlist.append([numch, startch, lastch])

#==========================
#== parallel writing pdf ==
#==========================
para_flag=1
# para_flag=0
#--
if para_flag==1:
    p=Pool(ncpus)
    list(p.map(write_pdf,inputlist))
    p.terminate()
else:
    list(map(write_pdf,inputlist))
