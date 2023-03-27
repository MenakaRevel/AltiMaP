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
        ix   = ix0 - 1
        iy   = iy0 - 1
        #---
        print (name, flag)
        #---------------------------
        lnames.append(name)
        l_lons.append(lon)
        l_lats.append(lat)
        lflags.append(flag)
#=============================      
lnames=np.array(lnames)
l_lons=np.array(l_lons)
l_lats=np.array(l_lats)
lflags=np.array(lflags)
# leledf=np.array(leledf)
# l_lons=np.array(l_lons)
# l_lats=np.array(l_lats)
N=float(len(lnames))
#=============================
mkdir("./fig")
mkdir("./fig/criteria")
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
# cmap=cm.rainbow_r
# #cmap.set_under("w",alpha=0)
# cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=1.0
vmax=8.0
norm=Normalize(vmin=vmin,vmax=vmax)

# bounds=np.arange(0.0,8.0,1.0)
# bounds=np.arange(0.0,6.0,1.0)
#-----------
# flag identity
# 10 = location was directly found
# 20 = location was on the unit-catchment outlet
# 30 = found the nearest river
# 31 = found the nearest main river
# 40 = correction for ocean grids
# 50 = bifurcation location
#cmap=colors.ListedColormap(['grey',"xkcd:ultramarine",'xkcd:clear blue','xkcd:jungle green',"xkcd:shamrock","xkcd:electric green","xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
#cmap=colors.ListedColormap(['grey','xkcd:jungle green',"xkcd:shamrock","xkcd:electric green","xkcd:ultramarine",'xkcd:clear blue',"xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
#cmap=colors.ListedColormap(['grey','xkcd:jungle green',"xkcd:shamrock","xkcd:greeny blue","xkcd:ultramarine",'xkcd:clear blue',"xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
# cmap=colors.ListedColormap(['grey',"xkcd:dark seafoam",'xkcd:deep teal',"xkcd:saffron","xkcd:purpleish",'xkcd:royal',"xkcd:peacock blue","xkcd:carmine","xkcd:black"])
# cmapL=matplotlib.colors.ListedColormap(['grey',"xkcd:dark seafoam",'xkcd:deep teal',"xkcd:saffron","xkcd:purpleish",'xkcd:royal',"xkcd:peacock blue","xkcd:carmine","xkcd:black"])
# marlist={10:'o', 20:'d', 30:'+', 31:'*', 32:'^', 40:'X', 50:'1'}
# corlist={10:'green', 20:'blue',30:'purple', 31:'yellow', 32:'xkcd:lavender', 40:'red',50:'xkcd:deep teal'}
# cmapL = matplotlib.colors.ListedColormap(['green', 'blue', 'purple', 'yellow','xkcd:lavender', 'red', 'xkcd:deep teal'])
marlist={10:'o', 20:'d', 30:'s',40:'^', 50:'P'}
corlist={10:'green', 20:'red', 30:'blue',  40:'yellow'} #,32:'xkcd:lavender',50:'xkcd:deep teal',30:'purple'}
cmapL = matplotlib.colors.ListedColormap(['green','red', 'blue', 'yellow'])

# cmapL.set_under("none") #"#000000",alpha=0)
# cmapL.set_over("none")
# cmapL.colorbar_extend="neither"
bounds=np.arange(0.0,5.0,1.0)
norml=BoundaryNorm(bounds,cmapL.N) #len(bounds)-1)

hgt= 11.69*(1.0/3.0)
wdt= 8.27
fig= plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(1,1)
ax = fig.add_subplot(G[0,0],projection=ccrs.Robinson())
#-----------------------------  
ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
#
## Draw the river network ##
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > tmp1.txt") 
# map(vec_par,np.arange(5,10+1,1))
map(vec_par,np.arange(2,10+1,1))
pnum=len(lnames)
for point in np.arange(pnum):
    # eled=leledf[point]
    name = lnames[point]
    lon  = l_lons[point]
    lat  = l_lats[point]
    flag = lflags[point]
    flag = flag_diff(flag)
    # c=cmapL(norml(flag))
    print (name,lon,lat,flag)
    c=corlist[flag]
    m=marlist[flag]
    ax.scatter(lon,lat,s=5,marker=m,zorder=110,edgecolors="k", facecolors=c, linewidth=0.1,transform=ccrs.PlateCarree()) #, 
#--
im=ax.scatter([],[],c=[],cmap=cmapL,s=0.1,vmin=vmin,vmax=vmax,norm=norml) # cmap=cmap, norm=norml
im.set_visible(False)
#cbar=M.colorbar(im,"right",size="2%")
ax.outline_patch.set_linewidth(0.0)
# colorbar
cax=fig.add_axes([0.40,0.10,0.35,.01])
cbar=plt.colorbar(im,orientation="horizontal",ticks=np.arange(0.5,8.0+1.0,1.0),cax=cax) #[10,20,30,31,40,50],extend='both',ticks=np.arange(0.0,1.0+0.001,0.1) extend='both',
cbar.set_ticklabels(['10', '20', '30', '40'])
# cbar.ax.set_ticklabels(['10', '20', '30', '31', '40', '50'])  # vertically oriented colorbar
cbar.ax.tick_params(labelsize=6)
cbar.set_label("Allocation Flags",fontsize=8)
plt.savefig("./fig/criteria/allocation_flag_map.png",dpi=500)
os.system("rm -r tmp*")