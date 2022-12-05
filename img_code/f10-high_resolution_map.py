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
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
tag="3sec"
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
# high-resolution data
visual=CaMa_dir+"/map/"+mapname+"/"+tag+"/s20w060.visual.bin"
visual=np.fromfile(visual,np.int8).reshape(12000,12000)
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

lllat = -20.
urlat = -19.75
lllon = -60.
urlon = -59.75

west  = lllon
east  = urlon
north = urlat
south = lllat

londiff=int((east-west)*1200)
latdiff=int((north-south)*1200)

npix= 0 #(90-north)*4
spix= latdiff #(90-south)*4
wpix= 0 #(180+west)*4
epix= londiff #(180+east)*4

#cmap=make_colormap(colors_list)
#cmap=mbar.colormap("H02")
# cmap=cm.seismic
# cmap=cm.rainbow_r
# #cmap.set_under("w",alpha=0)
# cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=1.0
vmax=26.0
norm=Normalize(vmin=vmin,vmax=vmax)

bounds=np.arange(-0.5,26,1.0)
#cmap=colors.ListedColormap(['grey',"xkcd:ultramarine",'xkcd:clear blue','xkcd:jungle green',"xkcd:shamrock","xkcd:electric green","xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
#cmap=colors.ListedColormap(['grey','xkcd:jungle green',"xkcd:shamrock","xkcd:electric green","xkcd:ultramarine",'xkcd:clear blue',"xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
#cmap=colors.ListedColormap(['grey','xkcd:jungle green',"xkcd:shamrock","xkcd:greeny blue","xkcd:ultramarine",'xkcd:clear blue',"xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
# cmap=colors.ListedColormap(['grey',"xkcd:dark seafoam",'xkcd:deep teal',"xkcd:saffron","xkcd:purpleish",'xkcd:royal',"xkcd:peacock blue","xkcd:carmine","xkcd:black"])
###
# sea=0
# land(undefined)=1
# land(defined in CaMa)=2 
# grid-box=3 
# catchment-boundary=5 
# channel=10 
# outlet-pixel=20 
# river-mouth=25
#
cmapL = matplotlib.colors.ListedColormap(['w','w','grey','k','w','k','w','y','w','w','blue','w','w','w','w','w','w','w','w','w','red', 'w','w','w','w','red'])
cmapL.set_under("none") #"#000000",alpha=0)
cmapL.set_over("none")
cmapL.colorbar_extend="neither"
norml=BoundaryNorm(bounds,cmapL.N) #len(bounds)-1)

# cmap=cmapL #cm.rainbow
# cmapL = cm.viridis_r
# norml=BoundaryNorm(bounds,cmapL.N)

print "making figure"
hgt= 11.69*(1.0/3.0)
wdt= 8.27
fig= plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(1,1)
ax = fig.add_subplot(G[0,0])#,projection=ccrs.Robinson())
#-----------------------------
m = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c',ax=ax)
#m.drawcoastlines( linewidth=0.1, color='k' )
# m.fillcontinents(color=land,lake_color=water,zorder=99)  
# ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
# ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
# #
m.drawparallels([lllat,urlat], labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
m.drawmeridians([lllon,urlon], labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
# pnum=len(lnames)
#--
print "imshow"
# [0:500,0:500]
im=m.imshow(ma.masked_less_equal(visual[npix:spix,wpix:epix],-9999.0),cmap=cmapL,norm=norml) # interpolation="nearest",origin="upper",
# im=ax.scatter([],[],c=[],cmap=cmapL,s=0.1,vmin=vmin,vmax=vmax,norm=norml) # cmap=cmap, norm=norml
# im.set_visible(False)
pnum=len(lnames)
for point in np.arange(pnum):
    lon  = l_lons[point]
    lat  = l_lats[point]
    ax.scatter(lon,lat,s=0.5,marker="o",zorder=110,edgecolors="g", facecolors="g")#,transform=ccrs.PlateCarree()) #, 
#-
#cbar=M.colorbar(im,"right",size="2%")
# ax.outline_patch.set_linewidth(0.0)
#colorbar
# cax=fig.add_axes([0.40,0.10,0.4,.01])
# cbar=plt.colorbar(im,orientation="horizontal",ticks=np.arange(0,25+0.1,1.0))#,cax=cax) #ticks=np.arange(vmin,vmax+0.1,1.0),,extend='both',ticks=np.arange(0.0,1.0+0.001,0.1) extend='both',
# cbar.ax.tick_params(labelsize=6)
# cbar.set_label("Visuals Flags",fontsize=8)
plt.savefig("./fig/criteria/high-resolution_map.png",dpi=500)