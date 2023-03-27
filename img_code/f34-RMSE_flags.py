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
eyear=2019
#=========================================
#TAG="CGLS"
TAG="HydroWeb"
# TAG="ICESat"
# TAG="HydroSat"
#=========================================
# odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
# fname = odir+"/hydroweb_cmf_daily_wse_VIC_BC.nc"
odir="/cluster/data6/menaka/Altimetry/results"
# odir = "/cluster/data6/menaka/CaMaVal/results"
if TAG=="HydroWeb":
    # fname=odir+"/HydroWeb/hydroweb_cmf_daily_wse_VIC_BC.nc"
    fname=odir+"/HydroWeb/hydroweb-hydroda_cmf_daily_wse_VIC_BC.nc"
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
print (pnum)
#print np.shape(sfcelv_hydroweb)
# colors=['xkcd:pastel blue','xkcd:aqua green','xkcd:soft pink']
labels=["$0km$","$<0.5km$","$<1km$","$<2km$","$<5km$","$<10km$","$>10km$"]
#==========
data0=[]
data1=[]
data2=[]
data3=[]
data4=[]
data5=[]
data6=[]
for point in np.arange(pnum):
    cmf_mean=np.mean(sfcelv_cmf[:,point])
    obs_mean=np.mean(ma.masked_equal(sfcelv_hydroweb[:,point],-9999.0))
    BIAS=abs(cmf_mean-obs_mean)
    RMSE1=RMSE(cmf_mean,obs_mean)
    # print (cmf_mean,obs_mean,BIAS)
    if disttomouth[point] <= 1.0e-5:
        data0.append(RMSE1)
    if disttomouth[point] <= 0.5:
        data1.append(RMSE1)
        print ("<=0.5km",pname[point], RMSE1, disttomouth[point])
    elif disttomouth[point] <= 1.0:
        data2.append(RMSE1)
        print ("<=1km",pname[point], RMSE1, disttomouth[point])
    elif disttomouth[point] <= 2.0:
        data3.append(RMSE1)
        print ("<=2km",pname[point], RMSE1, disttomouth[point])
    elif disttomouth[point] <= 5.0:
        data4.append(RMSE1)
        print ("<=5km",pname[point], RMSE1, disttomouth[point])
    elif disttomouth[point] <= 10.0:
        data5.append(RMSE1)
        print ("<=10km",pname[point], RMSE1, disttomouth[point])
    else:
        data6.append(RMSE1) 
        print (">10km",pname[point], RMSE1, disttomouth[point])
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
#--boxplot
#ax = fig.add_subplot(G[1,0])
ax = fig.add_subplot(G[0,0])
#ax.boxplot(data_area, labels=labels, showmeans=True)
flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
boxprops = dict(color='grey')#facecolor='none'
whiskerprops = dict(color='grey',linestyle="--")
capprops = dict(color='grey')
medianprops = dict(color='r')
meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
#with plt.style.context('classic'):#'default'):
print ("making box plot")
# print (np.mean(data,axis=0))
box=sns.boxplot(ax=ax,data=[data0,data1,data2,data3,data4,data5,data6]\
        ,fliersize=0.0, palette="husl", whis=1.5\
        ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
        ,meanprops=meanprops,capprops=capprops,medianprops=medianprops)
# box=ax.boxplot([data0,data1,data2,data3,data4,data5],labels=labels,meanline=True,patch_artist=True,showfliers=False)
# for patch, color in zip(box['boxes'], colors):
#         patch.set_facecolor(color)
# add number of Vss
for i,data in enumerate([data0,data1,data2,data3,data4,data5,data6]):
    nums="%d"%(len(data))
    print (labels[i],nums, pnum)
    # ax.text(i,np.median(data)+0.1*np.median(data),nums,ha="center",va="center",transform=ax.transAxes,fontsize=10)
ax.set_xticklabels(labels,rotation = 0)
ax.set_ylim(ymin=-0.2,ymax=13.2)
ax.set_ylabel('RMSE $(m)$', color='k',fontsize=8)
ax.set_xlabel('Distance to mouth $(km)$', color='k',fontsize=8)
ax.tick_params('y',labelsize=6, colors='k')
ax.tick_params('x',labelsize=6, colors='k')#,labelrotation=45)
ax.set_xticklabels(labels,rotation=0)
figname="RMSE_boxplot"
os.system("mkdir -p ./fig/"+TAG)
#--
print ("./fig/"+TAG+"/"+figname+".png")
plt.savefig("./fig/"+TAG+"/"+figname+".png",dpi=800,bbox_inches="tight", pad_inches=0.05)
