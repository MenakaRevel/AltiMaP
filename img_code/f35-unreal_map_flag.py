#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import datetime
import numpy as np
from numpy import ma
from numpy import dot
from numpy.linalg import solve
from numpy.polynomial.polynomial import Polynomial as P, polyvander as V
import re
import sys
import os
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
import warnings;warnings.filterwarnings('ignore')
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import errno

sys.path.append("./src")
from read_patchMS import upstream
from river_function import river_profile
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
#=============================
#=============================
colors=["xkcd:sea blue","xkcd:tangerine","xkcd:dark pink","xkcd:pinkish","xkcd:light urple"]
#======================================================================
rivername0="MEKONG" #"LENA" #"CONGO" #"AMAZONAS" #"AMAZONAS"
stream0="MEKONG" #"LENA" #"CONGO" #"AMAZONAS" #
# dataname="HydroWeb"
odir="/cluster/data6/menaka/Altimetry/fig/river_profile"
TAG="HydroWeb"
mapname="glb_06min"
restag="3sec"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
#======================================================================
mkdir("./fig")
mkdir("./fig/criteria")
#=============================
res=1.0/1200.0
nx =12000
ny =12000
hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
if restag == "3sec":
    res=1.0/1200.0
    nx =12000
    ny =12000
    hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
#======================================================================
# noobs="/cluster/data6/menaka/Altimetry/out/unreal_obs_20210622.txt"
noobs="/cluster/data6/menaka/Altimetry/out/unreal_obs_20210930.txt"
noVS=[]
with open(noobs,"r") as f:
    lines=f.readlines()
    for line in lines:
        line    = filter(None,re.split(" ",line))
        station = line[0].strip()
        noVS.append(station)
############################################################
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
kxlst=[]
kylst=[]
#===========================
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210709.txt"
obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210920.txt"
fname=obstxt
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
    eled    = float(line[7])
    EGM08   = float(line[8])
    EGM96   = float(line[9])
    sat     = line[10].strip()
    dist    = float(line[11])
    flag    = int(line[12])
    kx      = int(line[13])
    ky      = int(line[14])
    #-----------------------
    if station in noVS:
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
        kxlst.append(kx)
        kylst.append(ky)
#-----------------------------
pnum=len(nums)   
#=============================
# river width
sup=2
w=0.02
alpha=1
width=0.5

land="#C0C0C0"
water="#FFFFFF"

west=-180.0
east=180.0
north=90.0
south=-58.0

lllat = -58.
urlat = 90.
lllon = -180.
urlon = 180.

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

#cmap=make_colormap(colors_list)
#cmap=mbar.colormap("H02")
# cmap=cm.seismic
#cmap.set_under("w",alpha=0)
# cmap=cm.get_cmap("rainbow_r")
# cmapL=cmap #cm.get_cmap("rainbow_r")
# vmin=0.0
# vmax=2000.0
# norm=Normalize(vmin=vmin,vmax=vmax)
bounds=np.arange(0.0,8.0,1.0)
marlist={10:'o', 20:'d', 30:'+', 31:'*', 32:'^', 40:'X', 50:'1'}
corlist={10:'green', 20:'blue',30:'purple', 31:'yellow', 32:'xkcd:lavender', 40:'red',50:'xkcd:deep teal'}
cmapL = matplotlib.colors.ListedColormap(['green', 'blue', 'purple', 'yellow','xkcd:lavender', 'red', 'xkcd:deep teal'])
norml=BoundaryNorm(bounds,cmapL.N)
vmin=1.0
vmax=8.0
# norm=Normalize(vmin=vmin,vmax=vmax)

bounds=np.arange(0.0,8.0,1.0)
hgt=11.69*(1.0/3.0)
wdt=8.27
fig=plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(1,1)
ax=fig.add_subplot(G[0,0],projection=ccrs.Robinson())
#-----------------------------  
ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
#
pnum=len(nums)
for point in np.arange(pnum):
    flag=lflag[point]
    eled=leled[point]
    lon =lons[point]
    lat =lats[point]
    # c=cmapL(norml(flag))
    #print lon,lat,pname[point][0],mean_bias
    print (pname[point], flag)
    c=corlist[flag]
    m=marlist[flag]
    ax.scatter(lon,lat,s=0.5,marker="o",zorder=110,edgecolors=c, facecolors=c,transform=ccrs.PlateCarree())
#--
im=ax.scatter([],[],c=[],cmap=cmapL,s=0.1,vmin=vmin,vmax=vmax,norm=norml)#
im.set_visible(False)
#cbar=M.colorbar(im,"right",size="2%")
ax.outline_patch.set_linewidth(0.0)
#colorbar
cax=fig.add_axes([0.40,0.10,0.4,.01])
cbar=plt.colorbar(im,orientation="horizontal",ticks=np.arange(0.5,8.0+1.0,1.0),cax=cax) #,extend='both',extend='both',ticks=np.arange(0.0,1.0+0.001,0.1)
cbar.set_ticklabels(['10', '20', '30', '31', '32', '40', '50'])
cbar.ax.tick_params(labelsize=6)
cbar.set_label("Allocation Flag",fontsize=8)
plt.savefig("./fig/criteria/unrealstic_VS_map_flag.png",dpi=500)