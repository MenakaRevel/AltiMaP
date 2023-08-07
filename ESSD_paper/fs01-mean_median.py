#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cbook
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes,inset_axes,mark_inset
import matplotlib.patches as patches
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import re
import os
import errno
import pandas as pd
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
    try:
        with open(fname,"r") as f:
            lines=f.readlines()
    except:
        print("No file: "+fname)
        return -9999.0, -9999.0, -9999.0, -9999.0, -9999.0
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
        wse  = float(line[2]) +egm08-egm96
        data.append(wse)
    data=np.array(data)
    return np.mean(data), np.std(data), np.max(data), np.min(data), np.median(data)
#=============================
# sfcelv
syear=2002
eyear=2013
#=========================================
#TAG="CGLS"
TAG="HydroWeb"
# TAG="ICESat"
# TAG="HydroSat"
# # # #=========================================
# # # # odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
# # # # fname = odir+"/hydroweb_cmf_daily_wse_VIC_BC.nc"
# # # odir="/cluster/data6/menaka/Altimetry/results"
# # # if TAG=="HydroWeb":
# # #     fname0=odir+"/HydroWeb/hydroweb_cmf_daily_wse_VIC_BC.nc"
# # #     fname1=odir+"/HydroWeb/hydroweb_cmf_linear_daily_wse_VIC_BC.nc"
# # #     fname2=odir+"/HydroWeb/hydroweb_cmf_elediff_daily_wse_VIC_BC.nc"
# # #     fname3=odir+"/HydroWeb/hydroweb_cmf_updown_daily_wse_VIC_BC.nc"
# # # if TAG=="CGLS":
# # #     fname0=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
# # #     fname1=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
# # #     fname2=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
# # # if TAG=="ICESat":
# # #     fname0=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
# # #     fname1=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
# # #     fname3=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
# # # if TAG=="HydroSat":
# # #     fname0=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
# # #     fname1=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
# # #     fname2=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
# # # ############################################################
# # # ######################## original data #####################
# # # nc0 = xr.open_dataset(fname0)
# # # sfcelv_hydroweb0=nc0.sfcelv_hydroweb.values
# # # sfcelv_cmf0=nc0.sfcelv_cmf.values
# # # pname=nc0.name.values
# # # lons=nc0.lon.values
# # # lats=nc0.lat.values
# # # basins=nc0.Basin.values
# # # rivers=nc0.river.values
# # # countries=nc0.country.values
# # # sfcelv_hydroweb_max=nc0.sfcelv_hydroweb_max.values
# # # sfcelv_hydroweb_min=nc0.sfcelv_hydroweb_min.values
# # # sfcelv_cmf_max=nc0.sfcelv_cmf_max.values
# # # sfcelv_cmf_min=nc0.sfcelv_cmf_min.values
# # # disttomouth=nc0.disttomouth.values
# # # elediff=nc0.elediff.values
# # # flag=nc0.flag.values
# # # nc0.close()
# # # ######################## interpolated data #####################
# # # nc1 = xr.open_dataset(fname1)
# # # pname1=nc1.name.values
# # # sfcelv_cmf1=nc1.sfcelv_cmf.values
# # # nc1.close()
# # # ######################## elevation differnce data #####################
# # # nc2 = xr.open_dataset(fname2)
# # # pname2=nc2.name.values
# # # sfcelv_cmf2=nc2.sfcelv_cmf.values
# # # nc2.close()
# # # ######################## upstream downstream data #####################
# # # nc3 = xr.open_dataset(fname3)
# # # pname3=nc3.name.values
# # # sfcelv_cmf3=nc3.sfcelv_upstream_cmf.values
# # # sfcelv_cmf4=nc3.sfcelv_downstream_cmf.values
# # # nc3.close()
# # # ############################################################
mkdir("./fig")
# mkdir("./fig/"+TAG)
# mkdir("./fig/"+TAG+"/along_river_VS")
# pnum=10 #len(pname)
#print np.shape(sfcelv_hydroweb)
colors=['xkcd:pastel blue','xkcd:teal','xkcd:aqua green','xkcd:dark pink','xkcd:purple','xkcd:magenta']
labels=["cmf oroginal","cmf interpolated","cmf ele diff",TAG]
#=============================
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
tag="3sec"
res=1.0/1200.0
upthr=10.0
dwthr=10.0
thrs=15.0
#=============================
# Read the CMF variables
if mapname == 'glb_15min':
    nx      = 1440
    ny      = 720
    ny_     = 640
elif mapname == 'glb_06min':
    nx      = 3600
    ny      = 1800
    ny_     = 1500
elif mapname == 'glb_01min':
    nx      = 21600
    ny      = 10800
    ny_     = 10800
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
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
# rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
rivseq = np.fromfile(rivseq,np.int32).reshape(ny,nx)
#---
nextX=nextxy[0]
nextY=nextxy[1]
#=============================
# meanWSE_VICBC="/cluster/data6/menaka/ensemble_org/CaMa_out/GLBVIC_BC001/sfcelv_mean2000-2013.bin"
meanWSE_VICBC="/work/a06/menaka/ensemble_org/CaMa_out/GLBVIC_BC001/sfcelv_mean2000-2013.bin"
meanWSE_VICBC=np.fromfile(meanWSE_VICBC,np.float32).reshape(ny_,nx)
#=============================
vmin=1.0
vmax=26.0
norm=Normalize(vmin=vmin,vmax=vmax)
bounds=np.arange(-0.5,26,1.0)
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
#=============================
fname="../out/altimetry_"+mapname+"_20230407.txt"
df=pd.read_csv(fname,delim_whitespace=True)
print df.head()
#----------------
mean_elevation=[]
mean_WSE=[]
std_WSE=[]
range_WSE=[]
median_WSE=[]
cmf_WSE=[]
riv_hgt=[]
dist=0.0
# columns=df.columns
print df.columns
for point in df.index:
    meanW, stdW, maxW, minW, medianW = meanHydroWeb(df["station"][point],egm96=df["EGM96"][point],egm08=df["EGM08"][point])
    mean_WSE.append(meanW)
    std_WSE.append(stdW)
    range_WSE.append(maxW-minW)
    median_WSE.append(medianW)
#-------------------------------------------
hgt=11.69*(1.0/3.0)
wdt=8.27*(1.0/2.0)
fig=plt.figure(figsize=(wdt, hgt))
#plt.title(pname[point][0],fontsize=12)
G  = gridspec.GridSpec(1,1)
ax = fig.add_subplot(G[0,0])
ax.plot(mean_WSE,median_WSE,"o",color="grey",linestyle='none',linewidth=0.3,marker="o",fillstyle="none",markersize=2)
ax.plot([min(mean_WSE),max(mean_WSE)],[min(median_WSE),max(median_WSE)],"--",color="k",linestyle='solid',linewidth=0.2)
ax.set_xlim(0,max([max(mean_WSE),max(median_WSE)]))
ax.set_ylim(0,max([max(mean_WSE),max(median_WSE)]))
ax.set_xlabel("mean WSE [m]",fontsize=10)
ax.set_ylabel("median WSE [m]",fontsize=10)
plt.tight_layout()
plt.savefig("./fig/fs01-mean_vs_median.png",dpi=800)
plt.savefig("./fig/fs01-mean_vs_median.pdf",dpi=800)
plt.savefig("./fig/fs01-mean_vs_median.jpg",dpi=800)