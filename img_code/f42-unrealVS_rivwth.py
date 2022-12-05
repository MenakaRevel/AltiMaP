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
def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    where ever the observed data contains nan
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]

    return data[:,0],data[:,1]
#========================================
def RMSE(s,o):
    """
    Root Mean Squre Error
    input:
        s: simulated
        o: observed
    output:
        RMSE: Root Mean Squre Error
    """
    o=ma.masked_where(o<=0.0,o).filled(0.0)
    s=ma.masked_where(o<=0.0,s).filled(0.0)
    s,o = filter_nan(s,o)
    return np.sqrt(np.mean((s-o)**2))
#=========================================
def cname(lat,lon):
    if lat < 0.0:
        ns="s" 
    else:
        ns="n" 
    #-------
    if lon < 0.0:
        we="w" 
    else:
        we="e" 
    #############
    south="%02d"%(abs(int(math.floor(lat/10.0)*10)))
    west="%03d"%(abs(int(math.floor(lon/10.0)*10)))
    return ns+south+we+west
#=========================================
def westsouth(lat,lon):
    return float(int(math.floor(lon/10.0)*10)), float(int(math.floor(lat/10.0)*10))
#=========================================
def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5) / multiplier
#=============================
def round_half_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n*multiplier - 0.5) / multiplier
#=============================
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
# sfcelv
syear=2000
eyear=2020
start=datetime.date(syear,1,1)
end=datetime.date(eyear,12,31)
days=(end-start).days + 1
#=========================================
#TAG="CGLS"
TAG="HydroWeb"
# TAG="ICESat"
# TAG="HydroSat"
#=========================================
############################################################
# rivername0=sys.argv[1] #"CONGO" #"AMAZONAS"
dataname="HydroWeb" #sys.argv[2]
# odir=sys.argv[2] #"/cluster/data6/menaka/Altimetry/results" #sys.argv[3] #"/cluster/data6/menaka/Altimetry/results"
mapname="glb_06min" #sys.argv[4] #"glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4" #sys.argv[5] #"/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
restag="3sec" #sys.argv[6] #"3sec"
obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20220729.txt" #sys.argv[7] #"./out/altimetry_"+mapname+"_test.txt"
thr=10 #float(sys.argv[6])
stream0="AMAZONAS"
upthr = thr #10.0
dwthr = thr #10.0
#=============================
TAG=dataname
res=1.0/1200.0
nx =12000
ny =12000
hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
if restag == "3sec":
    res=1.0/1200.0
    nx =12000
    ny =12000
    hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
#=============================
# Read the CMF variables
if mapname == 'glb_15min':
    nXX     = 1440
    nYY     = 720
    nY_     = 640
    dXX     = 0
    dYY     = 0
elif mapname == 'glb_06min':
    nXX     = 3600
    nYY     = 1800
    nY_     = 1500
    dXX     = 0
    dYY     = 0
elif mapname == 'glb_01min':
    nXX     = 21600
    nYY     = 10800
    nY_     = 10800
    dXX     = 0
    dYY     = 0
elif mapname == 'amz_06min':
    nXX     = 350
    nYY     = 250
    nY_     = 250
    dXX     = 1000
    dYY     = 85
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
nextxy = np.fromfile(nextxy,np.int32).reshape(2,nYY,nXX)
# rivwth = np.fromfile(rivwth,np.float32).reshape(nYY,nXX)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(nYY,nXX)
# rivlen = np.fromfile(rivlen,np.float32).reshape(nYY,nXX)
elevtn = np.fromfile(elevtn,np.float32).reshape(nYY,nXX)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,nYY,nXX)
uparea = np.fromfile(uparea,np.float32).reshape(nYY,nXX)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(nYY,nXX)
rivseq = np.fromfile(rivseq,np.int32).reshape(nYY,nXX)
#---
nextX=nextxy[0]
nextY=nextxy[1]
#=============================
meanWSE_VICBC="/cluster/data6/menaka/ensemble_org/CaMa_out/GLBVIC_BC001/sfcelv_mean2000-2013.bin"
meanWSE_VICBC=np.fromfile(meanWSE_VICBC,np.float32).reshape(1500,3600)
############################################################
# pnum=10 #len(pname)
#print np.shape(sfcelv_hydroweb)
colors=['xkcd:pastel blue','xkcd:teal','xkcd:aqua green','xkcd:dark pink','xkcd:purple','xkcd:magenta']
labels=["cmf oroginal","cmf interpolated","cmf ele diff",TAG]
#=============================
vmin=1.0
vmax=26.0
norm=Normalize(vmin=vmin,vmax=vmax)
bounds=np.arange(-0.5,26,1.0)
############################################################
cmapL = matplotlib.colors.ListedColormap(['w','w','grey','k','w','k','w','y','w','w','blue','w','w','w','w','w','w','w','w','w','red', 'w','w','w','w','red'])
cmapL.set_under("none") #"#000000",alpha=0)
cmapL.set_over("none")
cmapL.colorbar_extend="neither"
norml=BoundaryNorm(bounds,cmapL.N) #len(bounds)-1)
############################################################
noobs="/cluster/data6/menaka/Altimetry/out/unreal_obs_20220729.txt"
noVS=[]
with open(noobs,"r") as f:
    lines=f.readlines()
    for line in lines:
        line    = filter(None,re.split(" ",line))
        station = line[0].strip()
        noVS.append(station)
#======================================================================
eldiff=[]
lrivwt=[]
urivwt=[]
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
fname=obstxt
#--
f=open(fname,"r")
lines=f.readlines()
for line in lines[1::]:
    line    = filter(None,re.split(" ",line))
    # print line
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
    rivwth  = float(line[19])
    #-----------------------
    # meanW, stdW = meanHydroWeb(station,egm96=EGM96,egm08=EGM08)
    # eldiff.append(meanW - elev)
    # west, south = westsouth(lat,lon)
    # north = south + 10.0
    # east  = west + 10.0
    # cname0 = cname(lat,lon)
    # # rivwidth
    # rivwth=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".rivwth.bin"
    # # print (rivwth)
    # rivwth=np.fromfile(rivwth,np.float32).reshape(12000,12000)
    lrivwt.append(rivwth)
    if station in noVS:
        urivwt.append(rivwth)
    # print (riv,station,kx,ky)
    # if riv != rivername0:
    #     continue
    # if stream != stream0:
    #     continue
    # nums.append(num)
    # river.append(riv)
    # pname.append(station)
    # lons.append(lon)
    # lats.append(lat)
    # xlist.append(ix)
    # ylist.append(iy)
    # lelev.append(elev)
    # egm08.append(EGM08)
    # egm96.append(EGM96)
    # llsat.append(sat)
    # ldtom.append(dist)
    # lflag.append(flag)
    # kxlst.append(kx)
    # kylst.append(ky)
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
sns.distplot(lrivwt,ax=ax, hist = True, kde = True,
    kde_kws = {'linewidth': 1,'linestyle':'-'},bins=50,
    label = "$River Width_{All}$", color="xkcd:coral",norm_hist=True)
sns.distplot(urivwt,ax=ax, hist = True, kde = True,
    kde_kws = {'linewidth': 1,'linestyle':'-'},bins=50,
    label = "$River Width_{unreal}$", color="xkcd:teal green",norm_hist=True)
# ax.set_xlim(xmin=-20.2,xmax=0.2)
ax.set_ylabel("density", color='k',fontsize=8)
ax.set_xlabel('River Width $(m)$', color='k',fontsize=8)
plt.savefig("./fig/criteria/unrealVS_rivwth.png",dpi=500)
print ("Mean: ",np.mean(lrivwt),"Median: ",np.median(lrivwt), "std: ",np.std(lrivwt) )    
print ("Mean: ",np.mean(urivwt),"Median: ",np.median(urivwt), "std: ",np.std(urivwt) )
