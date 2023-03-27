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
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
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
import seaborn as sns
import string

sys.path.append("../src")
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
#=================================
def vec_par(LEVEL):
    # river width
    ax=None
    sup=2
    w=0.005
    width=0.5
    ax=ax or plt.gca()
    txt=prename+"%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec "+prename+"00.txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    f = open(txt,"r")
    lines = f.readlines()
    f.close()
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        ix = int((lon1 + 180.)*4.0)
        iy = int((-lat1 + 90.)*4.0)

        if lon1-lon2 > 180.0:
            print lon1, lon2
            lon2=180.0
        elif lon2-lon1> 180.0:
            print lon1,lon2
            lon2=-180.0

        # if rivermap[iy,ix] == 1.0:
        colorVal="w" 
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#==============================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,transform=ccrs.PlateCarree(),linewidth=width,zorder=105,alpha=alpha)
#=============================
def flag_diff(flag):
    # flag identity
    # 10 = on the river centerline
    # 11 = on the river channel
    # 12 = location was on the unit-catchment outlet
    # 20 = found the nearest river
    # 21 = found the nearest main river
    # 30 = found the nearest perpendicular main river
    # 31 = bifurcation location
    # 40 = correction for ocean grids
    if flag == 10 or flag == 11 or flag == 12:
        return 10
    if flag == 20 or flag == 21:
        return 20
    if flag == 30 or flag == 31 or flag == 32:
        return 30
    if flag == 40:
        return 40
#=============================
def mk_boxplot(data,colors,ax=None):
    ax=ax or plt.gca()
    flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
    boxprops = dict(color='grey')#facecolor='none'
    whiskerprops = dict(color='grey',linestyle="--")
    capprops = dict(color='grey')
    medianprops = dict(color='grey',linestyle="-",linewidth=1.0)
    meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
    # make boplot
    box=sns.boxplot(ax=ax,data=[data], fliersize=0.0, palette=colors, whis=1.5\
        ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
        ,meanprops=meanprops,capprops=capprops,medianprops=medianprops) #"Paired"
    ax.set_xticklabels(["expert","ordinary"])
    # ax.set_xlabel(["expert","ordinary"])
    ax.set_ylim(ymin=-0.2,ymax=60.2)
    return 0
#=============================
def mk_hist(data,color="red",label="data",ax=None):
    ax=ax or plt.gca()
    sns.distplot(data,ax=ax, hist = True, kde = True,
            kde_kws = {'linewidth': 1,'linestyle':'-'},
            label = label, color=color, norm_hist=True)
    # ax.set_ylim(ymin=-0.2,ymax=60.2)
    ax.tick_params(axis='both', which='major', labelsize=6)
    ax.set_xlabel(label,fontsize=8)
    return 0 
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
prename=re.split("-",sys.argv[0])[0]
#=============================
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()
#=============================
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
# # Read the CMF variables
# if mapname == 'glb_15min':
#     nx      = 1440
#     ny      = 720
# elif mapname == 'glb_06min':
#     nx      = 3600
#     ny      = 1800
# elif mapname == 'glb_01min':
#     nx      = 21600
#     ny      = 10800
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
#======================================================================
mkdir("./fig")
# mkdir("./fig/criteria")
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
noobs="/cluster/data6/menaka/Altimetry/out/unreal_obs_20220729.txt"
noVS=[]
with open(noobs,"r") as f:
    lines=f.readlines()
    for line in lines:
        line    = filter(None,re.split(" ",line))
        station = line[0].strip()
        # print (station)
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
upare=[]
rivwth=[]
#===========================
obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20220730.txt"
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
    rw      = float(line[19])
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
        flag=flag_diff(flag)
        lflag.append(flag)
        kxlst.append(kx)
        kylst.append(ky)
        upare.append(uparea[iy,ix]*1e-6)
        rivwth.append(rw)
#-----------------------------
pnum=len(nums)   
#=============================
# river width
sup=2
w=0.01
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
bounds=np.arange(0.0,4.0,1.0)
# marlist={10:'o', 20:'d', 30:'+', 31:'*', 32:'^', 40:'X', 50:'1'}
# corlist={10:'green', 20:'blue',30:'purple', 31:'yellow', 32:'xkcd:lavender', 40:'red',50:'xkcd:deep teal'}
# cmapL = matplotlib.colors.ListedColormap(['green', 'blue', 'purple', 'yellow','xkcd:lavender', 'red', 'xkcd:deep teal'])

# marlist={10:'o', 20:'d', 30:'s',40:'^', 50:'P'}
# colors=["#4c72b0","#de8452","#55a868","#c54e51"]
# corlist={10:'#4c72b0', 20:'#de8452', 30:'#55a868', 40:'#c54e51'} #,32:'xkcd:lavender',50:'xkcd:deep teal',30:'purple'}
# cmapL = matplotlib.colors.ListedColormap(colors)

marlist={10:'o', 20:'d', 30:'s',40:'^'}
marsize={10:0.75, 20:0.75, 30:0.9,40:1.0}
colors=["#e0fbfb","#97c1d9","#3c5a80","#ef6c4d"]
corlist={10:'#e0fbfb', 20:'#97c1d9', 30:'#3c5a80', 40:'#ef6c4d'} #,32:'xkcd:lavender',50:'xkcd:deep teal',30:'purple'}
cmapL = matplotlib.colors.ListedColormap(colors)
zorder = {10:110, 20:111, 30:112, 40:115}
norml=BoundaryNorm(bounds,cmapL.N)
vmin=1.0
vmax=8.0
# norm=Normalize(vmin=vmin,vmax=vmax)

bounds=np.arange(0.0,8.0,1.0)
hgt=11.69*(1.0/2.0)
wdt=8.27
fig=plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(3,3)
ax1=fig.add_subplot(G[0:2,:],projection=ccrs.Robinson())
#-----------------------------  
ax1.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
#--
box="%f %f %f %f"%(lllon,urlon,urlat,lllat) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+"  > "+prename+"00.txt") 
# map(vec_par,np.arange(5,10+1,1))
map(vec_par,np.arange(2,10+1,1))
# map(vec_par,np.arange(5,10+1,1))
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
    ax1.scatter(lon,lat,s=5,marker=m,zorder=zorder[flag],edgecolors="k", facecolors=c, linewidth=0.15,transform=ccrs.PlateCarree())
#--
im=ax1.scatter([],[],c=[],cmap=cmapL,s=0.1,vmin=vmin,vmax=vmax,norm=norml)#
im.set_visible(False)
#cbar=M.colorbar(im,"right",size="2%")
ax1.outline_patch.set_linewidth(0.0)
# #colorbar
# cax=fig.add_axes([0.40,0.10,0.4,.01])
# cbar=plt.colorbar(im,orientation="horizontal",ticks=np.arange(0.5,8.0+1.0,1.0),cax=cax) #,extend='both',extend='both',ticks=np.arange(0.0,1.0+0.001,0.1)
# cbar.set_ticklabels(['10', '20', '30', '31', '32', '40', '50'])
# cbar.ax.tick_params(labelsize=6)
# cbar.set_label("Allocation Flag",fontsize=8)
ax1.text(0.00,1.00,"%s)"%(string.ascii_lowercase[0]),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
features=[]
for i,flag in enumerate([10,20,30,40]):
    features.append(mlines.Line2D([], [], color=corlist[flag], marker=marlist[flag],
                    markeredgecolor="k",markeredgewidth=0.5,markersize=2, label='%d'%(flag),linewidth=0.0))
l,b,w,h=0.55,0.40,0.45,0.01
# l,b,w,h=0.4,0.1,0.45,0.01
legend=plt.legend(handles=features,bbox_to_anchor=(l,b), loc="lower center",fontsize=8,
           bbox_transform=fig.transFigure, ncol=4,  borderaxespad=0.0, frameon=False)#
#====
ax2=fig.add_subplot(G[2,0])
mk_hist(np.log10(upare),color="red",label="$catchment$ $area$ $(km^2)$",ax=ax2)
print (np.median(upare))
ax2.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[1]),ha="left",va="center",transform=ax2.transAxes,fontsize=10)
ax2.set_xticks(np.arange(2,7+1,1))
ax2.set_xticklabels([r"$10^{%d}$"%(i) for i in np.arange(2,7+1,1)])
#====
ax3=fig.add_subplot(G[2,1])
mk_hist(leled,color="green",label="$elevation$ $(m)$",ax=ax3)
print (np.median(leled))
ax3.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[2]),ha="left",va="center",transform=ax3.transAxes,fontsize=10)
#====
ax4=fig.add_subplot(G[2,2])
mk_hist(rivwth,color="xkcd:slate",label="$river$ $width$ $(m)$",ax=ax4)
print (np.median(rivwth))
ax4.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[3]),ha="left",va="center",transform=ax4.transAxes,fontsize=10)
#=====
plt.savefig("./fig/f05-unrealstic_VS_map_flag.png",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./fig/f05-unrealstic_VS_map_flag.jpg",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./fig/f05-unrealstic_VS_map_flag.pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
os.system("rm -r "+prename+"*.txt")