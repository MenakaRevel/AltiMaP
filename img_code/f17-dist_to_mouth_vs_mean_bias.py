#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
from numpy.lib.function_base import median
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
from scipy import stats
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
# hwb_amptd=ma.masked_where(sfcelv_hydroweb_max==-9999.0,sfcelv_hydroweb_max-sfcelv_hydroweb_min).filled(-9999.0)
# cmf_amptd=ma.masked_where(sfcelv_cmf_max==-9999.0,sfcelv_cmf_max-sfcelv_cmf_min).filled(-9999.0)
pnum=len(pname)
#print np.shape(sfcelv_hydroweb)
colors=['xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
labels=["$<1m$","$<5m$","$>5m$"]
#==========
# masked out -9999 values
# sfcelv_cmf=ma.masked_where(sfcelv_hydroweb==-9999.0,sfcelv_cmf).filled(-9999.0)
sfcelv_bias=ma.masked_where(sfcelv_hydroweb==-9999.0,(sfcelv_cmf-sfcelv_hydroweb)).filled(-9999.0)#.compressed()#
mean_bias=np.mean(ma.masked_where(sfcelv_bias==-9999.0,sfcelv_bias),axis=0)
mean_bias=np.abs(mean_bias)
mean_bias=ma.masked_greater(mean_bias,50.0)
#----
for point in np.arange(pnum):
    mean_bias0=np.mean(ma.masked_where(sfcelv_bias[:,point]==-9999.0,sfcelv_bias[:,point]))
    lon=lons[point]
    lat=lats[point]
    if np.abs(mean_bias0) > 10.0:
        print (lon,lat,pname[point][0],mean_bias0)
#----
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
# scatter
ax0 = fig.add_subplot(G[0,0])
ax0.plot(disttomouth,mean_bias,color="k",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=2)
#-- fit line --
slope, intercept, r_value, p_value, std_err = stats.linregress(disttomouth,np.abs(mean_bias))
print (slope, intercept, r_value, p_value, std_err)
ax0.plot(disttomouth,slope*disttomouth+intercept,color="r",linestyle='-',linewidth=1)
ax0.set_title("All")
ax0.set_ylabel('Long-term Bias $(m)$', color='k',fontsize=8)
ax0.set_xlabel('Distance to mouth $(km)$', color='k',fontsize=8)
ax0.set_xlim(xmin=0.0)
ax0.set_ylim(ymin=0.0)
ax0.tick_params('y',labelsize=6, colors='k')
ax0.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
#----
figname="Bias_vs_disttomouth_scatter"
os.system("mkdir -p ./fig/"+TAG)
#--
print ("./fig/"+TAG+"/"+figname+".png")
plt.savefig("./fig/"+TAG+"/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)