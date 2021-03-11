#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
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
from matplotlib.colors import LogNorm,Normalize,ListedColormap
import matplotlib.cm as cm
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature

sys.path.append("./src")
# import read_hydroweb as hweb
import read_cgls as cgls
import read_hydrosat as hsat
import read_icesat as isat
import read_grrats as grt
from read_upstream import upstream
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
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# print sys.argv
# mapname=sys.argv[1]
# CaMa_dir=sys.argv[2]
# ncpus=int(sys.argv[3])
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
#---
nextX=nextxy[0]
nextY=nextxy[1]
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
syear=2003
eyear=2020
start=datetime.date(syear,1,1)
end=datetime.date(eyear,12,31)
days=(end-start).days + 1
#=============================
lnames=[]
leledf=[]
l_lons=[]
l_lats=[]
lupare=[]
lelevt=[]
#=============================
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
        ix0  = int(line[5]) 
        iy0  = int(line[6])
        eled = float(line[7])
        egm08= float(line[8])
        egm96= float(line[9])
        sat  = line[10].strip()
        flag = int(line[11])
        ix   = ix0 - 1
        iy   = iy0 - 1
        #---------------------------
        if flag == 1:
            lnames.append(name)
            leledf.append(abs(eled))
            l_lons.append(lon)
            l_lats.append(lat)
            lupare.append(uparea[iy,ix]*1e-6)
            lelevt.append(elevtn[iy,ix])
        else:
            continue

#=============================      
lnames=np.array(lnames)
leledf=np.array(leledf)
l_lons=np.array(l_lons)
l_lats=np.array(l_lats)
lupare=np.array(lupare)
lelevt=np.array(lelevt)
N=float(len(lnames))
#=============================
mkdir("./fig")
mkdir("./fig/criteria")
#=============================
hgt=11.69*(1.0/3.0)
wdt=8.27
fig=plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(1,1)
ax = fig.add_subplot(G[0,0])
#-----------------------------    
ax.plot(leledf,lupare,marker="o",fillstyle="none",markerfacecolor=None,markeredgecolor="red",markersize=3,linestyle=None,linewidth=0.0,alpha=1)
# ax.set_xscale('log')
# scentific notaion
# ax.ticklabel_format(style="sci",axis="x",scilimits=(0,0))
# ax.xaxis.major.formatter._useMathText=True 
# ax.set_xticks(flags)
# ax.set_xticklabels(flags,fontsize=8)
ax.set_ylabel('Upstream Catchment Area $(km^2)$', color='k',fontsize=8)
ax.set_xlabel('Elevation Difference', color='k',fontsize=8)
ax.tick_params('y',labelsize=6, colors='k')
ax.tick_params('x',labelsize=6, colors='k')
plt.savefig("./fig/criteria/ele_diff_uparea.png",dpi=500)
#=============================
plt.close()
hgt=11.69*(1.0/3.0)
wdt=8.27
fig=plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(1,1)
ax = fig.add_subplot(G[0,0])
#-----------------------------    
ax.plot(leledf,lelevt,marker="o",fillstyle="none",markerfacecolor=None,markeredgecolor="red",markersize=3,linestyle=None,linewidth=0.0,alpha=1)
# ax.set_xscale('log')
# scentific notaion
# ax.ticklabel_format(style="sci",axis="x",scilimits=(0,0))
# ax.xaxis.major.formatter._useMathText=True 
# ax.set_xticks(flags)
# ax.set_xticklabels(flags,fontsize=8)
ax.set_ylabel('Elevation $(m)$', color='k',fontsize=8)
ax.set_xlabel('Elevation Difference', color='k',fontsize=8)
ax.tick_params('y',labelsize=6, colors='k')
ax.tick_params('x',labelsize=6, colors='k')
plt.savefig("./fig/criteria/ele_diff_elevtn.png",dpi=500)