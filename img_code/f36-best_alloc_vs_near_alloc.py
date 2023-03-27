#!/opt/local/bin/python
# -*- coding: utf-8 -*-
'''
Compare the best allocation vs nearst allocation
'''
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
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
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
def flag_diff(flag):
    if flag == 10 or flag == 20:
        return 10
    if flag == 30 or flag == 32:
        return 20
    if flag == 31 or flag == 50:
        return 30
    if flag == 40:
        return 40
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
        data.append(wse)
    data=np.array(data)
    return np.mean(data), np.std(data)
#=====================================
#=============================
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# print sys.argv
# mapname=sys.argv[1]
# CaMa_dir=sys.argv[2]
# ncpus=int(sys.argv[3])
#=============================
mkdir("./fig")
mkdir("./fig/criteria")
dataname="HydroWeb"
odir="/cluster/data6/menaka/Altimetry/results"
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
restag="3sec"
res=1.0/1200.0
nx =12000
ny =12000
hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
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
kx1lt=[]
ky1lt=[]
kx2lt=[]
ky2lt=[]
elev1=[]
elev2=[]
muwse=[]
#-------------------------------------------
# fname="./out/altimetry_"+mapname+"_test.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210618.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210701.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210706.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210909.txt"
fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210920.txt"
############################################################
f=open(fname,"r")
lines=f.readlines()
inputlist=[]
index=0
for line in lines[1::]:
    line    = filter(None,re.split(" ",line))
    #print line
    num     = line[0]
    station = line[1].strip()
    line2   = re.split("_",station)
    # print num, line2
    riv     = line2[1]
    stream  = line2[2]
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
    kx1     = int(line[13])
    ky1     = int(line[14])
    kx2     = int(line[15])
    ky2     = int(line[16])
    #---
    # if flag != 31: #or flag != 50:
    if flag not in [31,32,50]:
        continue
    if kx2==-9999.0 or ky2==-9999.0:
        continue
    if kx1==kx2 and ky1==ky2:
        continue
    inputlist.append([str(index),station,str(lon),str(lat),str(EGM08),str(EGM96),str(kx1),str(ky1),str(kx2),str(ky2)])
    index=index+1
    pname.append(station)
    lons.append(lon)
    lats.append(lat)
    egm08.append(EGM08)
    egm96.append(EGM96)
    kx1lt.append(kx1)
    ky1lt.append(ky1)
    kx2lt.append(kx2)
    ky2lt.append(ky2)
    print (riv, stream, flag)
pnum=len(pname)

# multiprocessing array
elev1=np.ctypeslib.as_ctypes(np.zeros([pnum],np.float32))
shared_array_elev1  = sharedctypes.RawArray(elev1._type_, elev1)
elev2=np.ctypeslib.as_ctypes(np.zeros([pnum],np.float32))
shared_array_elev2  = sharedctypes.RawArray(elev2._type_, elev2)
muwse=np.ctypeslib.as_ctypes(np.zeros([pnum],np.float32))
shared_array_muwse  = sharedctypes.RawArray(muwse._type_, muwse)

def read_data(inputlist):
    index=int(inputlist[0])
    station=inputlist[1]
    lon=float(inputlist[2])
    lat=float(inputlist[3])
    EGM08=float(inputlist[4])
    EGM96=float(inputlist[5])
    kx1=int(inputlist[6])
    ky1=int(inputlist[7])
    kx2=int(inputlist[8])
    ky2=int(inputlist[9])
    #=====
    tmp_elev1  = np.ctypeslib.as_array(shared_array_elev1)
    tmp_elev2  = np.ctypeslib.as_array(shared_array_elev2)
    tmp_muwse  = np.ctypeslib.as_array(shared_array_muwse)
    #=====
    cname0 = cname(lat,lon)
    elevtn0 = CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".elevtn.bin"
    # print (visual)
    elevtn0 = np.fromfile(elevtn0,np.float32).reshape(12000,12000)
    #===
    tmp_elev1[index]=elevtn0[ky1-1,kx1-1]
    tmp_elev2[index]=elevtn0[ky2-1,kx2-1]
    # read HydroWeb data
    meanW, stdW = meanHydroWeb(station,egm96=EGM96,egm08=EGM08)
    tmp_muwse[index]=meanW
    print (riv, stream, meanW, elevtn0[ky1-1,kx1-1], elevtn0[ky2-1,kx2-1])
###
p     = Pool(20)
res   = p.map(read_data, inputlist)
elev1 = np.ctypeslib.as_array(shared_array_elev1)
elev2 = np.ctypeslib.as_array(shared_array_elev2)
muwse = np.ctypeslib.as_array(shared_array_muwse)
p.terminate()
############################################################
hgt= 11.69*(1.0/2.0)
wdt= 8.27
fig= plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(2,2)
ax1 = fig.add_subplot(G[0,0])
ax1.plot(elev1,muwse,color="xkcd:sea blue",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
ax1.set_ylabel('Mean WSE $(m)$', color='k',fontsize=8)
ax1.set_xlabel('Best Allocation: Elevation  $(m)$', color='k',fontsize=8)
ax1.tick_params('y',labelsize=6, colors='k')
ax1.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
ax1.plot([0.0,4500.0],[0.0,4500.0],linestyle="--",linewidth=0.5,color="k")
ax1.set_xlim(xmin=0.0,xmax=4500.0)
ax1.set_ylim(ymin=0.0,ymax=4500.0)
# inset axes....
axins1 = zoomed_inset_axes(ax1, zoom=2, loc='lower right')
# fix the number of ticks on the inset axes
axins1.yaxis.get_major_locator().set_params(nbins=10)
axins1.xaxis.get_major_locator().set_params(nbins=10)

plt.setp(axins1.get_xticklabels(), visible=False)
plt.setp(axins1.get_yticklabels(), visible=False)

# axins1 = ax1.inset_axes([0.0, 0.5, 0.47, 0.47])
axins1.plot(elev1,muwse,color="xkcd:sea blue",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
axins1.plot([0.0,1000.0],[0.0,1000.0],linestyle="--",linewidth=0.5,color="k")
# sub region of the original image
x1, x2, y1, y2 = 0.0, 1000.0, 0.0, 1000.0
axins1.set_xlim(x1, x2)
axins1.set_ylim(y1, y2)
axins1.set_xticklabels('')
axins1.set_yticklabels('')

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax1, axins1, loc1=2, loc2=4, fc="none", ec="0.5")

#============================
ax2 = fig.add_subplot(G[0,1])
ax2.plot(elev2,muwse,color="xkcd:dark turquoise",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
# ax2.set_ylabel('Mean WSE $(m)$', color='k',fontsize=8)
ax2.set_xlabel('Nearest Allocation: Elevation  $(m)$', color='k',fontsize=8)
ax2.tick_params('y',labelsize=6, colors='k')
ax2.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
ax2.plot([0.0,4500.0],[0.0,4500.0],linestyle="--",linewidth=0.5,color="k")
ax2.set_xlim(xmin=0.0,xmax=4500.0)
ax2.set_ylim(ymin=0.0,ymax=4500.0)
# inset axes....
axins2 = zoomed_inset_axes(ax2, zoom=2, loc='lower right')
# fix the number of ticks on the inset axes
axins2.yaxis.get_major_locator().set_params(nbins=10)
axins2.xaxis.get_major_locator().set_params(nbins=10)

plt.setp(axins2.get_xticklabels(), visible=False)
plt.setp(axins2.get_yticklabels(), visible=False)

# axins2 = ax2.inset_axes([0.0, 0.5, 0.47, 0.47])
axins2.plot(elev1,muwse,color="xkcd:dark turquoise",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
axins2.plot([0.0,1000.0],[0.0,1000.0],linestyle="--",linewidth=0.5,color="k")
# sub region of the original image
x1, x2, y1, y2 = 0.0, 1000.0, 0.0, 1000.0
axins2.set_xlim(x1, x2)
axins2.set_ylim(y1, y2)
axins2.set_xticklabels('')
axins2.set_yticklabels('')

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax2, axins2, loc1=2, loc2=4, fc="none", ec="0.5")

#============================
ax3 = fig.add_subplot(G[1,0])
ax3.plot(elev1,elev2,color="xkcd:reddish pink",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
ax3.set_ylabel('Nearest Allocation: Elevation  $(m)$', color='k',fontsize=8)
ax3.set_xlabel('Best Allocation: Elevation  $(m)$', color='k',fontsize=8)
ax3.tick_params('y',labelsize=6, colors='k')
ax3.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
ax3.plot([0.0,4500.0],[0.0,4500.0],linestyle="--",linewidth=0.5,color="k")
ax3.set_xlim(xmin=0.0,xmax=4500.0)
ax3.set_ylim(ymin=0.0,ymax=4500.0)
# inset axes....
axins3 = zoomed_inset_axes(ax3, zoom=2, loc='lower right')
# fix the number of ticks on the inset axes
axins3.yaxis.get_major_locator().set_params(nbins=10)
axins3.xaxis.get_major_locator().set_params(nbins=10)

plt.setp(axins3.get_xticklabels(), visible=False)
plt.setp(axins3.get_yticklabels(), visible=False)
# axins3 = ax3.inset_axes([0.0, 0.5, 0.47, 0.47])
axins3.plot(elev1,muwse,color="xkcd:reddish pink",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
axins3.plot([0.0,1000.0],[0.0,1000.0],linestyle="--",linewidth=0.5,color="k")
# sub region of the original image
x1, x2, y1, y2 = 0.0, 1000.0, 0.0, 1000.0
axins3.set_xlim(x1, x2)
axins3.set_ylim(y1, y2)
axins3.set_xticklabels('')
axins3.set_yticklabels('')

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax3, axins3, loc1=2, loc2=4, fc="none", ec="0.5")

plt.savefig("./fig/criteria/scatter_best_nearset.png",dpi=500)
plt.savefig("./fig/criteria/scatter_best_nearset.pdf",dpi=500)