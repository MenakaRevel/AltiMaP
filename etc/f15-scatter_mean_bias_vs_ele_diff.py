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
dataxx0=[]
datayy0=[]
dataxx1=[]
datayy1=[]
dataxx2=[]
datayy2=[]
dataxx3=[]
datayy3=[]
for point in np.arange(pnum):
    cmf_mean=np.mean(sfcelv_cmf[:,point])
    obs_mean=np.mean(ma.masked_equal(sfcelv_hydroweb[:,point],-9999.0))
    BIAS=abs(cmf_mean-obs_mean)
    # print (cmf_mean,obs_mean,BIAS)
    datayy0.append(BIAS)
    dataxx0.append(abs(elediff[point]))
    if disttomouth[point] <= 1.0:
        datayy1.append(BIAS)
        dataxx1.append(abs(elediff[point]))
        print ("<=1m",pname[point], BIAS, disttomouth[point])
    elif disttomouth[point] <= 5.0:
        datayy2.append(BIAS)
        dataxx2.append(abs(elediff[point]))
        print ("<=5m",pname[point], BIAS, disttomouth[point])
    else:
        datayy3.append(BIAS)
        dataxx3.append(abs(elediff[point]))
        print (">5m",pname[point], BIAS, disttomouth[point])
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
G = gridspec.GridSpec(2,2)

# scatter
ax0 = fig.add_subplot(G[0,0])
ax0.plot(dataxx0,datayy0,color="k",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
ax0.set_title("All")

ax1 = fig.add_subplot(G[0,1])
ax1.plot(dataxx0,datayy0,color="k",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
ax1.set_title("<1m")

ax2 = fig.add_subplot(G[1,0])
ax2.plot(dataxx1,datayy1,color="k",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
ax2.set_title("<5m")

ax3 = fig.add_subplot(G[1,1])
ax3.plot(dataxx2,datayy2,color="k",linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)
ax3.set_title(">5m")
#----
figname="Bias_scatter"
os.system("mkdir -p ./fig/"+TAG)
#--
print ("./fig/"+TAG+"/"+figname+".png")
plt.savefig("./fig/"+TAG+"/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
