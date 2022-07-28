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
from matplotlib.colors import LogNorm,Normalize,ListedColormap
import matplotlib.cm as cm
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
# colors=['xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
# labels=["$<1km$","$<5km$","$>5km$"]
#==========
data0=[]
for point in np.arange(pnum):
    cmf_mean=np.mean(sfcelv_cmf[:,point])
    obs_mean=np.mean(ma.masked_equal(sfcelv_hydroweb[:,point],-9999.0))
    BIAS=abs(cmf_mean-obs_mean)
    RMSE1=RMSE(cmf_mean,obs_mean)
    data0.append(RMSE1)
#----
#cmap=make_colormap(colors_list)
#cmap=mbar.colormap("H02")
cmap=cm.viridis_r
#cmap.set_under("w",alpha=0)
cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=0.0
vmax=10.0
norm=Normalize(vmin=vmin,vmax=vmax)
#-----------
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
G = gridspec.GridSpec(1,1)
#ax = fig.add_subplot(G[1,0])
ax = fig.add_subplot(G[0,0])
for point in np.arange(pnum):
    dist=np.abs(disttomouth[point])
    elev=np.abs(elediff[point])
    c=cmapL(norm(data0[point]))
    #print lon,lat,pname[point][0],mean_bias
    # ax.scatter(dist,elev,s=0.5,marker="o",zorder=110,edgecolors=c, facecolors=none)
    ax.plot(dist,elev,color=c,linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
im=ax.scatter([],[],c=[],cmap=cmapL,s=0.1,vmin=vmin,vmax=vmax,norm=norm)#
im.set_visible(False)
ax.set_ylabel('Elevation differnce $(m)$', color='k',fontsize=8)
ax.set_xlabel('Distance to mouth $(km)$', color='k',fontsize=8)
ax.tick_params('y',labelsize=6, colors='k')
ax.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
# ax.set_xticklabels(labels,rotation=0)
#colorbar
# cax=fig.add_axes([0.40,0.10,0.4,.01])
cbar=plt.colorbar(im) #,orientation="horizontal",extend='max',ticks=np.arange(vmin,vmax+0.1,1.0),cax=cax) #,extend='both',ticks=np.arange(0.0,1.0+0.001,0.1)
cbar.ax.tick_params(labelsize=6)
cbar.set_label("RMSE $(m)$",fontsize=8)

figname="RMSE_scatter_disttom_elediff"
os.system("mkdir -p ./fig/"+TAG)
#--
print ("./fig/"+TAG+"/"+figname+".png")
plt.savefig("./fig/"+TAG+"/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)