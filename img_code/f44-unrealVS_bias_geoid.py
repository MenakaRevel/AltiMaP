#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
from numpy.core.fromnumeric import mean
matplotlib.use('Agg')
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import re
import sys
import os
import errno
from scipy import stats
from sklearn.metrics import r2_score

sys.path.append("./src")
from read_patchMS import upstream
from river_function import river_profile
import read_hydroweb as hweb
import read_cgls as cgls
import read_hydrosat as hsat
import read_icesat as isat
import read_grrats as grt
#=========================================
def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
#========================================
def meanHydroWeb(station,egm08=0.0,egm96=0.0): 
    fname="/cluster/data6/menaka/HydroWeb/data/hydroprd_"+station+".txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    head=33
    #--
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
        wse  = float(line[2]) #+egm08-egm96
        if wse >= 9999.0:
            continue
        data.append(wse)
    data=np.array(data)
    return np.mean(data), np.std(data)
#=====================================

############################################################
# rivername0=sys.argv[1] #"CONGO" #"AMAZONAS"
dataname="HydroWeb" #sys.argv[2]
# odir=sys.argv[2] #"/cluster/data6/menaka/Altimetry/results" #sys.argv[3] #"/cluster/data6/menaka/Altimetry/results"
mapname="glb_06min" #sys.argv[4] #"glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4" #sys.argv[5] #"/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
restag="3sec" #sys.argv[6] #"3sec"
obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20220729.txt" #sys.argv[7] #"./out/altimetry_"+mapname+"_test.txt"
thr=10 #float(sys.argv[6])
stream0="AMAZONAS"
upthr = thr #10.0
dwthr = thr #10.0
############################################################
noobs="/cluster/data6/menaka/AltiMaP/out/unreal_obs_20220802.txt"
noVS=[]
bias=[]
with open(noobs,"r") as f:
    lines=f.readlines()
    for line in lines:
        line    = filter(None,re.split(" ",line))
        station = line[0].strip()
        bias_   = float(line[6].strip())
        noVS.append(station)
        bias.append(bias_)
#======================================================================
############################################################
nums=[]
river=[]
pname=[]
lons =[]
lats =[]
xlist=[]
ylist=[]
lelev=[]
egm08=[]
egm96=[]
llsat=[]
ldtom=[]
lflag=[]
kxlst=[]
kylst=[]
#-------------------------------------------
# fname="./out/altimetry_"+mapname+"_test.txt"
# fname="./out/altimetry_"+mapname+"_20210518.txt"
obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20220725.txt"
fname=obstxt
#--
f=open(fname,"r")
lines=f.readlines()
for line in lines[1::]:
    line    = filter(None,re.split(" ",line))
    #print line
    num     = line[0]
    station = line[1]
    line2   = re.split("_",station)
    riv     = line2[1]
    stream  = line2[2]
    lon     = float(line[3])
    lat     = float(line[4])
    ix      = int(line[5])-1
    iy      = int(line[6])-1
    elev    = float(line[7])
    EGM08   = float(line[8])
    EGM96   = float(line[9])
    sat     = line[10].strip()
    dist    = float(line[11])
    flag    = int(line[12])
    kx      = int(line[13])
    ky      = int(line[14])
    # rivwth  = float(line[19])
    #-----------------------
    # print (riv,station,kx,ky)
    # if riv != rivername0:
    #     continue
    # if stream != stream0:
    #     continue
    # if station not in noVS:
    #     continue
    nums.append(num)
    river.append(riv)
    pname.append(station)
    lons.append(lon)
    lats.append(lat)
    xlist.append(ix)
    ylist.append(iy)
    lelev.append(elev)
    egm08.append(EGM08)
    egm96.append(EGM96)
    llsat.append(sat)
    ldtom.append(dist)
    lflag.append(flag)
    kxlst.append(kx)
    kylst.append(ky)
    # print (riv,station)
#=============================
pnum=len(pname)
#=====================================
ldiff=[]
for point in np.arange(pnum):
    meanW, stdW = meanHydroWeb(pname[point],egm96=egm96[point],egm08=egm08[point])
    iXX = xlist[point]
    iYY = ylist[point]
    elev= lelev[point]
    sat = llsat[point]
    diff= meanW - elev
    ldiff.append(diff)
    if abs(egm08[point]) > 10.0 and abs(diff) > 40.0:
        if abs(diff) > 1.2*abs(egm08[point]):
            print (pname[point], diff, egm08[point])
#=============================
mkdir("./fig")
mkdir("./fig/criteria")
#=============================
hgt=11.69*(1.0/3.0)
wdt=8.27
fig=plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(1,1)
ax = fig.add_subplot(G[0,0])
#=============================
ax.plot(np.abs(egm08),np.abs(ldiff),marker="o",markersize=2.0, markerfacecolor="None",markeredgecolor='b',markeredgewidth=0.5,linewidth=0)
ax.plot([0.0,100.0],[0.0,100.0],linewidth=1,linestyle="--",color="grey")
ax.set_xlim(xmin=0.0,xmax=100.0)
ax.set_ylim(ymin=0.0,ymax=100.0)
ax.set_ylabel("Elevation difference [MERIT - Mean WSE] $(m)$", color='k',fontsize=8)
ax.set_xlabel('Geoid height [EGM08] $(m)$', color='k',fontsize=8)
plt.savefig("./fig/criteria/unrealVS_geoid.png",dpi=500)