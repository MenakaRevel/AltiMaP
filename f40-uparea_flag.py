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
import pandas as pd
import matplotlib.patches as mpatches
import re
import os
import sys
import errno
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from matplotlib import colors
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
def flag_diff(flag):
    if flag == 10 or flag == 20:
        return 10
    if flag == 30 or flag == 32:
        return 20
    if flag == 31 or flag == 50:
        return 30
    if flag == 40:
        return 40
#==================================
def vec_par(LEVEL,ax=None,sup=2,w=0.005,width=0.5):
    '''
    Plotting the river network
    '''
    ax=ax or plt.gca()
    txt="tmp_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec tmp1.txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print LEVEL, width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    with open(txt,"r") as f:
        lines = f.readlines()

    #print LEVEL, width, lines, txt
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        ix = int((lon1 - west)*(1/gsize))
        iy = int((-lat1 + north)*(1/gsize))

        if rivermap[iy,ix] == 0:
            continue

        if lon1-lon2 > 180.0:
            # print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            # print (lon1,lon2)
            lon2=-180.0
        #--------
        colorVal="w" #"grey"#
        #print (lon1,lon2,lat1,lat2,width)
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#==================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
#==================================
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# print sys.argv
# mapname=sys.argv[1]
# CaMa_dir=sys.argv[2]
# ncpus=int(sys.argv[3])
#----
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
lon0   = float(filter(None, re.split(" ",lines[4]))[0])
lat0   = float(filter(None, re.split(" ",lines[7]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])
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
rivermap=(nextxy[0]>0)*1.0
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
lflags=[]
# leledf=[]
l_lons=[]
l_lats=[]
#=============================
# fname="./out/altimetry_"+mapname+"_test.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210618.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210701.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210706.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210909.txt"
fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210920.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_glb_06min_20210920.txt"
upas={}
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
        flag = int(line[12])
        flag = flag_diff(flag)
        ix   = ix0 - 1
        iy   = iy0 - 1
        #---
        print (name, flag)
        if flag not in upas.keys():
            upas[flag]=[uparea[iy,ix]*1e-6]
        else:
            upas[flag].append(uparea[iy,ix]*1e-6)
#----
# colors=["xkcd:reddy brown","xkcd:dark pink","xkcd:pinkish","xkcd:light urple","xkcd:tangerine"]
colors=["#4c72b0","#de8452","#55a868","#c54e51"]
# plt.clf()  
# plt.close() 
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(1.0/2.0)
fig=plt.figure(figsize=(wdt,hgt))
G = gridspec.GridSpec(1,1)
ax = fig.add_subplot(G[0,0])
#=====
# df=pd.DataFrame(upas.items(), columns=['10','20','30','40'])
# df=pd.DataFrame.from_dict(upas)
# print (df.head())
# print type(upas[10]),type(upas[20]),type(upas[30]),type(upas[40])
# sns.kdeplot(ax=ax,data=df,cumulative=True, common_norm=False, common_grid=True) ##
# box=sns.boxplot(ax=ax,data=[upas[10],upas[20],upas[30],upas[40]])
# ax.set_xticklabels(labels,rotation = 0)
# ax.set_ylim(ymin=-0.2,ymax=13.2)
labels=["Flag 10","Flag 20","Flag 30","Flag 40"]
for i in np.arange(0,4,1):
    sns.distplot(np.log(upas[(i+1)*10]),ax=ax, hist = False, kde = True,
            kde_kws = {'linewidth': 1,'linestyle':'-'},
            label = labels[i],color=colors[i],norm_hist=True) 
ax.set_ylabel('density', color='k',fontsize=8)
ax.set_xlabel('$log(upstream catchment area)$ $(km)$', color='k',fontsize=8)
ax.tick_params('y',labelsize=6, colors='k')
ax.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
plt.savefig("./fig/criteria/uparea_flag.png",dpi=500)
plt.savefig("./fig/criteria/uparea_flag.pdf",dpi=500)
plt.savefig("./fig/criteria/uparea_flag.jpg",dpi=500)