#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
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
#=========================================
#=============================
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
# sfcelv
syear=2002
eyear=2013
#=========================================
#TAG="CGLS"
TAG="HydroWeb"
# TAG="ICESat"
# TAG="HydroSat"
#=========================================
# odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
# fname = odir+"/hydroweb_cmf_daily_wse_VIC_BC.nc"
odir="/cluster/data6/menaka/Altimetry/results"
if TAG=="HydroWeb":
    fname0=odir+"/HydroWeb/hydroweb_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/HydroWeb/hydroweb_cmf_linear_daily_wse_VIC_BC.nc"
    fname2=odir+"/HydroWeb/hydroweb_cmf_elediff_daily_wse_VIC_BC.nc"
if TAG=="CGLS":
    fname0=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
    fname2=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
if TAG=="ICESat":
    fname0=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
    fname3=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
if TAG=="HydroSat":
    fname0=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
    fname2=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
############################################################
######################## original data #####################
nc0 = xr.open_dataset(fname0)
sfcelv_hydroweb0=nc0.sfcelv_hydroweb.values
sfcelv_cmf0=nc0.sfcelv_cmf.values
pname=nc0.name.values
lons=nc0.lon.values
lats=nc0.lat.values
basins=nc0.Basin.values
rivers=nc0.river.values
countries=nc0.country.values
sfcelv_hydroweb_max=nc0.sfcelv_hydroweb_max.values
sfcelv_hydroweb_min=nc0.sfcelv_hydroweb_min.values
sfcelv_cmf_max=nc0.sfcelv_cmf_max.values
sfcelv_cmf_min=nc0.sfcelv_cmf_min.values
disttomouth=nc0.disttomouth.values
elediff=nc0.elediff.values
flag=nc0.flag.values
nc0.close()
######################## interpolated data #####################
nc1 = xr.open_dataset(fname1)
pname1=nc1.name.values
sfcelv_cmf1=nc1.sfcelv_cmf.values
nc1.close()
######################## elevation differnce data #####################
nc2 = xr.open_dataset(fname2)
pname2=nc2.name.values
sfcelv_cmf2=nc2.sfcelv_cmf.values
nc2.close()
############################################################
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
mkdir("./fig")
mkdir("./fig/"+TAG)
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
vmax=4.0
norm=Normalize(vmin=vmin,vmax=vmax)

bounds=np.arange(0.5,4.0,1.0)
#cmap=colors.ListedColormap(['grey',"xkcd:ultramarine",'xkcd:clear blue','xkcd:jungle green',"xkcd:shamrock","xkcd:electric green","xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
#cmap=colors.ListedColormap(['grey','xkcd:jungle green',"xkcd:shamrock","xkcd:electric green","xkcd:ultramarine",'xkcd:clear blue',"xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
#cmap=colors.ListedColormap(['grey','xkcd:jungle green',"xkcd:shamrock","xkcd:greeny blue","xkcd:ultramarine",'xkcd:clear blue',"xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
# cmap=colors.ListedColormap(['grey',"xkcd:dark seafoam",'xkcd:deep teal',"xkcd:saffron","xkcd:purpleish",'xkcd:royal',"xkcd:peacock blue","xkcd:carmine","xkcd:black"])
cmapL = matplotlib.colors.ListedColormap(['green', 'blue','red']) #'purple', 'yellow', ])
# cmapL=matplotlib.colors.ListedColormap(['grey',"xkcd:dark seafoam",'xkcd:deep teal',"xkcd:saffron","xkcd:purpleish",'xkcd:royal',"xkcd:peacock blue","xkcd:carmine","xkcd:black"])
cmapL.set_under("none") #"#000000",alpha=0)
cmapL.set_over("none")
cmapL.colorbar_extend="neither"
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
pnum=len(pname)
#=============================
for point in np.arange(pnum):
    #prepare
    cmf0=sfcelv_cmf0[:,point]
    cmf1=sfcelv_cmf1[:,point]
    cmf2=sfcelv_cmf2[:,point]
    org=sfcelv_hydroweb0[:,point]
    #----------------
    RMSE0=RMSE(cmf0,org)
    RMSE1=RMSE(cmf1,org)
    RMSE2=RMSE(cmf2,org)
    minRMSE=np.argmin(np.array([RMSE0,RMSE1,RMSE2])) + 1
    lon  = lons[point]
    lat  = lats[point]
    c=cmapL(norml(minRMSE))
    print (lon,lat,minRMSE)
    ax.scatter(lon,lat,s=0.5,marker="o",zorder=110,edgecolors=c, facecolors=c,transform=ccrs.PlateCarree()) #, 
#--
im=ax.scatter([],[],c=[],cmap=cmapL,s=0.1,vmin=vmin,vmax=vmax,norm=norml) # cmap=cmap, norm=norml
im.set_visible(False)
#cbar=M.colorbar(im,"right",size="2%")
ax.outline_patch.set_linewidth(0.0)
#colorbar
cax=fig.add_axes([0.40,0.10,0.4,.01])
cbar=plt.colorbar(im,orientation="horizontal",ticks=np.arange(vmin,vmax+0.1,1.0),cax=cax) #,extend='both',ticks=np.arange(0.0,1.0+0.001,0.1) extend='both',
cbar.ax.tick_params(labelsize=6)
cbar.set_label("RMSE comaprison",fontsize=8)
plt.savefig("./fig/"+TAG+"/RMSE_compare_map.png",dpi=500)