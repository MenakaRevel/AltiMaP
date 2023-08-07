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
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
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
from matplotlib import colors
import matplotlib.cm as cm
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
import pandas as pd
import warnings;warnings.filterwarnings('ignore')

sys.path.append("../src")

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
#==================================
def vec_par(LEVEL):
    '''
    Plotting the river network
    '''
    ax=None;sup=2;w=0.005;width=0.5
    ax=ax or plt.gca()
    txt=prename+"%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec "+prename+"00.txt 1 "+str(LEVEL)+" > "+txt)
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
        colorVal="w" #"grey" # "k"  #
        # print (lon1,lon2,lat1,lat2,width)
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#==================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    alpha=1.0
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,transform=ccrs.PlateCarree(),zorder=105,alpha=alpha)
#==================================
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
prename=re.split("-",sys.argv[0])[0]
# print sys.argv
# mapname=sys.argv[1]
# CaMa_dir=sys.argv[2]
# ncpus=int(sys.argv[3])
#=============================
# Read the CMF variables
if mapname == 'glb_15min':
    nx      = 1440
    ny      = 720
    gsize   = 0.25
elif mapname == 'glb_06min':
    nx      = 3600
    ny      = 1800
    gsize   = 0.1
elif mapname == 'glb_01min':
    nx      = 21600
    ny      = 10800
    gsize   = 1/60.
#=============================
fname="/work/a07/uddin/Cama-Flood/CaMa-Flood_v4/map/GBM_15min/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()
#-------
nxx    = int(filter(None, re.split(" ",lines[0]))[0])
nyy    = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
lon0   = float(filter(None, re.split(" ",lines[4]))[0])
lat0   = float(filter(None, re.split(" ",lines[7]))[0])
west   = float(filter(None, re.split(" ",lines[4]))[0])
east   = float(filter(None, re.split(" ",lines[5]))[0])
south  = float(filter(None, re.split(" ",lines[6]))[0])
north  = float(filter(None, re.split(" ",lines[7]))[0])
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
# nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
# rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
# rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
# rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
# elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
# lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
# nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
# rivseq = np.fromfile(rivseq,np.int32).reshape(ny,nx)
#---
# nextX=nextxy[0]
# nextY=nextxy[1]
# #---
# rivermap=(nextX>0)*1.0
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
# fname="./out/altimetry_"+mapname+"_test.txt"
# fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210618.txt"
# fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210701.txt"
# fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210706.txt"
# fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210909.txt"
# fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210920.txt"
# fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20220725.txt" # Used for ESSD
# fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20230407.txt"
fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20230727.txt" # CryoSat2
############################################################
df=pd.read_csv(fname, sep='\s+', header=0) #,encoding=str)
print (df.head())
print (df.columns)
df["uparea"]=[uparea[iy-1,ix-1]*1e-6 for ix, iy in zip(df["ix"],df["iy"])]
N=float(len(df))
#=============================
flag_fig={}
fflag=[]
fsta=[]
flon=[]
flat=[]
fkx1=[]
fky1=[]
fkx2=[]
fky2=[]
# fpox=[0,0,0,0,1,2,3]
# fpoy=[0,1,2,3,3,3,3]
fpox=[0,1,2,3]
fpoy=[3,3,3,3]
#=============================
mkdir("./fig")
# mkdir("./fig/criteria")
dataname="HydroWeb"
odir="/cluster/data6/menaka/AltiMaP/results"
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
restag="3sec"
res=1.0/1200.0
nXX =12000
nYY =12000
hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
#=============================
rivnum="/cluster/data6/menaka/CaMa_subbasin/output/rivnum_"+mapname+".bin"
print ("rivnum",rivnum)
rivnum=np.fromfile(rivnum,np.int32).reshape(ny,nx)
rivermap=((nextxy[0]>0)*(rivnum<=1))*1.0
#=============================
# river width
sup=2
w=0.02
alpha=1
width=0.5
#=============================
# color defintion for lan and ocean
land="#C0C0C0"
water="#FFFFFF"

# land="k"
land="w"

# west=-180.0
# east=180.0
# north=90.0
# south=-58.0
north=min(north,31.5)
lllat = north
urlat = south
lllon = west
urlon = east

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

#cmap=make_colormap(colors_list)
#cmap=mbar.colormap("H02")
# cmap=cm.seismic
# cmap=cm.rainbow_r
# #cmap.set_under("w",alpha=0)
# cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=1.0
vmax=8.0
norm=Normalize(vmin=vmin,vmax=vmax)

# bounds=np.arange(0.0,8.0,1.0)
# bounds=np.arange(0.0,6.0,1.0)
bounds=np.arange(0.0,4.0,1.0)
#-----------
# flag identity
# 10 = on the river centerline
# 11 = on the river channel
# 12 = location was on the unit-catchment outlet
# 20 = found the nearest river
# 21 = found the nearest main river
# 30 = found the nearest perpendicular main river
# 31 = bifurcation location
# 40 = correction for ocean grids
#======================================
marlist={10:'o', 20:'d', 30:'s',40:'^'}
marsize={10:3.5, 20:3.5, 30:5.5,40:7.5}
# colors=["#fcb17c","#dfc5fe","#218833","#fb020e"]
# corlist={10:'#fcb17c', 20:'#dfc5fe', 30:'#218833', 40:'#fb020e'} #,32:'xkcd:lavender',50:'xkcd:deep teal',30:'purple'}
colors=["#e0fbfb","#97c1d9","#3c5a80","#ef6c4d"]
corlist={10:'#e0fbfb', 20:'#97c1d9', 30:'#3c5a80', 40:'#ef6c4d'} #,32:'xkcd:lavender',50:'xkcd:deep teal',30:'purple'}
# New colors based on chat-gpt
colors = [
    "#CD5C5C",  # Terracotta
    "#228B22",  # Forest Green
    "#FFD700",  # Sunflower Yellow
    "#4169E1"  # Royal Blue
]
corlist={10:'#CD5C5C', 20:'#228B22', 30:'#FFD700', 40:'#4169E1'}
cmapL = matplotlib.colors.ListedColormap(colors)
zorder = {10:110, 20:111, 30:112, 40:115}
# cmapL.set_under("none") #"#000000",alpha=0)
# cmapL.set_over("none")
# cmapL.colorbar_extend="neither"
norml=BoundaryNorm(bounds,cmapL.N) #len(bounds)-1)
############################################################
vmin=1.0
vmax=26.0
normM=Normalize(vmin=1,vmax=26)
boundsM=np.arange(-0.5,26,1.0)
cmapM = matplotlib.colors.ListedColormap(['w','w','grey','k','w','k','w','y','w','w','blue','w','w','w','w','w','w','w','w','w','red', 'w','w','w','w','red'])
cmapM.set_under("none") #"#000000",alpha=0)
cmapM.set_over("none")
cmapM.colorbar_extend="neither"
normm=BoundaryNorm(boundsM,cmapM.N) #len(bounds)-1)
############################################################
hgt= 11.69*(2.0/5.0)
wdt= 8.27
fig= plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(1,1)
ax = fig.add_subplot(G[0,0],projection=ccrs.Robinson())
#-----------------------------  
ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', edgecolor='face', facecolor=land),zorder=100)
#
ax.add_geometries(Reader("/work/a06/menaka/Uddin_ERL/GBM_basin/GBM_basin.shp").geometries(),
                  ccrs.PlateCarree(),
                  facecolor='k', edgecolor='black',zorder=100)
## Draw the river network ##
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > "+prename+"00.txt") 
# map(vec_par,np.arange(8,10+1,1))
# map(vec_par,np.arange(5,10+1,1))
map(vec_par,np.arange(2,10+1,1))
#
pnum=len(df)

# df["corlist"]=[corlist[math.floor(flag/10.0)*10.0] for flag in df["flag"]]
# df["marlist"]=[marlist[math.floor(flag/10.0)*10.0] for flag in df["flag"]]
# df["marsize"]=[marsize[math.floor(flag/10.0)*10.0] for flag in df["flag"]]
# df["zorder"]=[zorder[math.floor(flag/10.0)*10.0] for flag in df["flag"]]

df["corlist"]=[corlist[math.floor((flag-100)/10.0)*10.0 if flag >= 100 else math.floor(flag/10.0)*10.0] for flag in df["flag"]]
df["marlist"]=[marlist[math.floor((flag-100)/10.0)*10.0 if flag >= 100 else math.floor(flag/10.0)*10.0] for flag in df["flag"]]
df["marsize"]=[marsize[math.floor((flag-100)/10.0)*10.0 if flag >= 100 else math.floor(flag/10.0)*10.0] for flag in df["flag"]]
df["zorder"]=[zorder[math.floor((flag-100)/10.0)*10.0 if flag >= 100 else math.floor(flag/10.0)*10.0] for flag in df["flag"]]

print (df.head())
# df.plot(x="lon", y="lat", s="marsize", marker=df["marlist"].tolist(), zorder=df["zorder"], edgecolors="k", facecolors=df["corlist"], linewidth=0.05, transform=ccrs.PlateCarree())
# # df.plot(x="lon", y="lat", s="marsize", marker=df["marlist"], zorder=df["zorder"], edgecolors="k", facecolors=df["corlist"], linewidth=0.05, transform=ccrs.PlateCarree())
for point in np.arange(len(df)):
    ax.scatter(df["lon"][point],df["lat"][point],s=df["marsize"][point],marker=df["marlist"][point],zorder=df["zorder"][point],edgecolors="k", facecolors=df["corlist"][point],linewidth=0.05,transform=ccrs.PlateCarree()) #,
#--
# imL=ax.scatter([],[],c=[],cmap=cmapL,s=0.1,vmin=vmin,vmax=vmax,norm=norml) # cmap=cmap, norm=norml
# imL.set_visible(False)
#cbar=M.colorbar(im,"right",size="2%")
ax.outline_patch.set_linewidth(0.0)

#================
# # add a historam
# l,b,w,h=ax.get_position().bounds
# ax_hist=fig.add_axes([l-0.1,b+0.05,w*(1.0/4.0),h*(1.0/3.0)])
# labels=["Flag 10","Flag 20","Flag 30","Flag 40"]
# for i in np.arange(0,4,1):
#     # print (labels[i], np.median(upas[(i+1)*10]))
#     sns.distplot(np.log10(df["uparea"][df["flag"]==(i+1)*10]),ax=ax_hist, hist = False, kde = True,
#             kde_kws = {'linewidth': 1,'linestyle':'-'},
#             label = labels[i],color=colors[i],norm_hist=True) 
# ax_hist.set_ylabel('density', color='k',fontsize=8)
# ax_hist.set_xlabel('$upstream$ $catchment$ $area$ $(km^2)$', color='k',fontsize=8)
# ax_hist.tick_params('y',labelsize=6, colors='k')
# ax_hist.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
# ax_hist.set_xticks(np.arange(1,9+1,1))
# ax_hist.set_xticklabels([r"$10^{%d}$"%(i) for i in np.arange(1,9+1,1)])
# ax_hist.spines['top'].set_visible(False)
# ax_hist.spines['right'].set_visible(False)
# ax_hist.set_facecolor(color="None")
#================
# legend
features=[]
for i,flag in enumerate([10,20,30,40]):
    # features.append(mlines.Line2D([], [], color=corlist[flag], marker=marlist[flag],
    #                 markeredgecolor="k",markeredgewidth=0.5,markersize=4,label='%d'%(flag),linewidth=0.0))
    features.append(mpatches.Patch(color=corlist[flag],label='%d'%(flag)))
l,b,w,h=0.4,0.1,0.45,0.01
legend=plt.legend(handles=features,bbox_to_anchor=(l,b), loc="lower left",
           bbox_transform=fig.transFigure, ncol=4,  borderaxespad=0.0, frameon=False)#

plt.tight_layout()

plt.savefig("./fig/fs02-allocation_flag_map_CryoSat2.png",dpi=800)
plt.savefig("./fig/fs02-allocation_flag_map_CryoSat2.pdf",dpi=800)
plt.savefig("./fig/fs02-allocation_flag_map_CryoSat2.jpg",dpi=800)

os.system("rm -r "+prename+"*.txt")