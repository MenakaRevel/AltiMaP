#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import matplotlib.patches as mpatches
import os

os.system("ln -sf ./src/read_patchMS.so read_patchMS.so")
from read_patchMS import upstream
#=========================================
# functions
#=========================================
def slope(ix,iy,nextxy,uparea,elevtn,rivseq):
    if rivseq[iy,ix]>1:
        nextX=nextxy[0]
        nextY=nextxy[1]
        uXX, uYY = upstream(ix+1,iy+1,nextX.T,nextY.T,uparea.T)
        uXX = uXX - 1
        uYY = uYY - 1
        slp=(elevtn[uYY,uXX]-elevtn[iy,ix])
    else:
        slp=0
    return slp
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
#=========================================
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
#=========================================
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
    fname=odir+"/HydroWeb/hydroweb_cmf_daily_wse_VIC_BC.nc"
if TAG=="CGLS":
    fname=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
if TAG=="ICESat":
    fname=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
if TAG=="HydroSat":
    fname=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
nc = xr.open_dataset(fname)
sfcelv_hydroweb=nc.sfcelv_hydroweb.values
sfcelv_cmf=nc.sfcelv_cmf.values
pname=nc.name.values
lons=nc.lon.values
lats=nc.lat.values
xvals=nc.variables["x-coord"].values
yvals=nc.variables["y-coord"].values
basins=nc.Basin.values
rivers=nc.river.values
countries=nc.country.values
sfcelv_hydroweb_max=nc.sfcelv_hydroweb_max.values
sfcelv_hydroweb_min=nc.sfcelv_hydroweb_min.values
sfcelv_cmf_max=nc.sfcelv_cmf_max.values
sfcelv_cmf_min=nc.sfcelv_cmf_min.values
disttomouth=nc.disttomouth.values
elediff=nc.elediff.values
nc.close()
pnum=len(pname)
#print np.shape(sfcelv_hydroweb)
colors=['xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
labels=["$<1km$","$<5km$","$>5km$"]
#==========
data0=[]
data1=[]
data2=[]
data3=[]
for point in np.arange(pnum):
    cmf_mean=np.mean(sfcelv_cmf[:,point])
    obs_mean=np.mean(ma.masked_equal(sfcelv_hydroweb[:,point],-9999.0))
    BIAS=abs(cmf_mean-obs_mean)
    ix = xvals[point]
    iy = yvals[point]
    elevtion=elevtn[iy,ix]
    # slp=slope(ix,iy,nextxy,uparea,elevtn,rivseq)
    slp=abs(elediff[point])
    RMSE1=RMSE(sfcelv_cmf[:,point],sfcelv_hydroweb[:,point])
    ## slope > 10 m --> dist_to_mouth <  1 km
    ## slope >  5 m --> dist_to_mouth <  5 km
    ## slope >  1 m --> dist_to_mouth <  10 km
    ## slope <  1 m --> dist_to_mouth <  All
    # print (cmf_mean,obs_mean,BIAS)
    if slp > 10.0 and disttomouth[point] <= 1.0:
        data0.append(RMSE1)
        print ("slope: ",slp,"<=1km", disttomouth[point],pname[point], RMSE1)
    elif slp > 5.0 and disttomouth[point] <= 5.0:
        data1.append(RMSE1)
        print ("slope: ",slp,"<=5km", disttomouth[point],pname[point], RMSE1)
    elif slp > 1.0 and disttomouth[point] <= 10.0:
        data2.append(RMSE1)
        print ("slope: ",slp,"<=5km", disttomouth[point],pname[point], RMSE1)
    elif slp <= 1.0 and disttomouth[point] < 1e20:
        data3.append(RMSE1)
        print ("slope: ",slp,"<=5km", disttomouth[point],pname[point], RMSE1)
    else: 
        print ("Not belong to criteria --->","elevation: ",slp,"distance", disttomouth[point],pname[point][0], RMSE1)
#----
colors=['xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink','xkcd:tan']
labels=["elevation$>10m$\ndistance$<1km$","elevation$>5m$\ndistance$<5km$","elevation$>1m$\ndistance$>10km$","elevation$<1m$\ndistance$:All$"]
# plt.clf()  
# plt.close() 
# figure in A4 size
va_margin= 0.0#1.38#inch 
ho_margin= 0.0#1.18#inch
hgt=(11.69 - 2*va_margin)*(1.0/3.0)
wdt=(8.27 - 2*ho_margin)*(1.0/2.0)
#fig=plt.figure(figsize=(8.27,11.69))
fig=plt.figure(figsize=(wdt,hgt))
#fig.suptitle("auto-correalated area")#"maximum autocorrelation length")
#G = gridspec.GridSpec(2,1)
G = gridspec.GridSpec(2,2)
#--boxplot
#ax = fig.add_subplot(G[1,0])
#ax.boxplot(data_area, labels=labels, showmeans=True)
flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
boxprops = dict(color='grey')#facecolor='none'
whiskerprops = dict(color='grey',linestyle="--")
capprops = dict(color='grey')
medianprops = dict(color='r')
meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
#with plt.style.context('classic'):#'default'):
print ("making box plot")
ax = fig.add_subplot(G[0,0])
# print (np.mean(data,axis=0))
box=ax.boxplot([data0,data1,data2,data3],labels=labels,meanline=True,patch_artist=True,showfliers=False)
for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
ax.set_ylabel('RMSE $(m)$', color='k',fontsize=8)
ax.set_xlabel('Distance to mouth $(km)$', color='k',fontsize=8)
ax.tick_params('y',labelsize=6, colors='k')
ax.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
ax.set_xticklabels(labels,fontsize=4,rotation=90)
#-----
figname="Bias_boxplot_criteria_elevation"
os.system("mkdir -p ./fig/"+TAG)
#--
print ("./fig/"+TAG+"/"+figname+".png")
plt.savefig("./fig/"+TAG+"/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)