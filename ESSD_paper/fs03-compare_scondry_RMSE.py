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
import xarray as xr
import math
import seaborn as sns
import re
import os
import sys
import errno
import string
import pandas as pd
from multiprocessing import Pool
from multiprocessing import Process
import warnings;warnings.filterwarnings('ignore')

sys.path.append("../src")
from read_sfcelv import read_sfcelv, read_sfcelv_multi
import read_hydroweb as hweb
#==============================
# statistic functions
#==============================
def filter_nan(s,o):
    """
    this functions removed the data  from simulated and observed data
    where ever the observed data contains nan
    """
    data = np.array([s.flatten(),o.flatten()])
    data = np.transpose(data)
    data = data[~np.isnan(data).any(1)]

    return data[:,0],data[:,1]
#==============================
def rmse(s,o):
    s,o = filter_nan(s,o)
    s   = ma.masked_where(o<=-9000.0,s).compressed() #filled(-9999.0)
    o   = ma.masked_where(o<=-9000.0,o).compressed() #filled(-9999.0)
    return np.sqrt(np.mean(ma.masked_where(o<=-9000.0,(s-o)**2)))
#==============================
def bias(s,o):
    s,o = filter_nan(s,o)
    s   = ma.masked_where(o<=-9000.0,s).compressed() #filled(-9999.0)
    o   = ma.masked_where(o<=-9000.0,o).compressed() #filled(-9999.0)
    return np.mean(ma.masked_where(o<=-9000.0,s))-np.mean(ma.masked_where(o<=-9000.0,o))
#=================================
def corr(s,o):
    s,o = filter_nan(s,o)
    s   = ma.masked_where(o<=-9000.0,s).compressed()
    o   = ma.masked_where(o<=-9000.0,o).compressed()
    return np.corrcoef(s,o)[0,1]
#======================================================================================================
# functions to read wse
#======================================================================================================
def read_wse_multi(ix, iy, syear, eyear, number, lat, lon):
    #print ix1,iy1
    wse = np.zeros( (len(ix), nbdays), 'f')
    wse_max = np.zeros( (len(ix), nbyears), 'f')
    wse_min = np.zeros( (len(ix), nbyears), 'f')
    wse_max_loc = np.zeros( (len(ix), nbyears), 'f')
    wse_min_loc = np.zeros( (len(ix), nbyears), 'f')
    for year in range(syear, eyear+1):
        #print year
        s_days = int( (datetime.date(year , 1,1) - datetime.date(syear, 1, 1)). days)
        e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)

        f = runoff_folder + '/sfcelv'+str(year)+'.bin'
        #print f, len(ix), len(iy)
        tmp = read_sfcelv_multi( ix, iy, e_days-s_days, f, nx, ny)

        #print year, e_days - s_days, s_days, e_days, outflw.shape
        wse[:,s_days:e_days] = tmp
        wse_max[:,year-syear] = np.nanmax(tmp, axis=1)
        wse_min[:,year-syear] = np.nanmin(tmp, axis=1)
        wse_max_loc[:,year-syear] = np.argmax(tmp, axis=1)
        wse_min_loc[:,year-syear] = np.argmin(tmp, axis=1)

    return wse , wse_max, wse_min, wse_max_loc, wse_min_loc
#======================================================================================================
# functions to read wse
#======================================================================================================
# read WSE from CMF outputs
def read_wse(ix, iy, syear, eyear, number, lat, lon):
    wse = np.zeros( nbdays, 'f')
    wse_max = np.zeros( nbyears, 'f')
    wse_min = np.zeros( nbyears, 'f')
    wse_max_loc = np.zeros( nbyears, 'f')
    wse_min_loc = np.zeros( nbyears, 'f')
    for year in range(syear, eyear+1):
        s_days = int( (datetime.date(year , 1,1) - datetime.date(syear, 1, 1)). days)
        e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)

        f = runoff_folder + '/sfcelv'+str(year)+'.bin'

        #outflw = np.fromfile(f, 'float32').reshape(-1,len(lat_global), len(lon_global))
        #if ix2 < 0 and iy2 < 0:
        #    tmp = outflw[:,iy1, ix1]
        #else:
        #    tmp = outflw[:,iy1, ix1] + outflw[:,iy2, ix2]

        #print tmp
        #sno=len(numbers)
        tmp = read_sfcelv( ix, iy, e_days-s_days, f, nx, ny)

        #print year, e_days - s_days, s_days, e_days, outflw.shape
        wse[s_days:e_days] = tmp
        wse_max[year-syear] = np.nanmax(tmp)
        wse_min[year-syear] = np.nanmin(tmp)
        wse_max_loc[year-syear] = np.argmax(tmp)
        wse_min_loc[year-syear] = np.argmin(tmp)

        '''
        # this is to check if the location is right or not.
        fig = plt.figure(figsize=(ssize ,ssize*0.5),dpi=300)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        plt.title('GRDC: '+ str(number))

        tmp = np.nanmean(outflw, axis=0)
        tmp = np.ma.masked_where(tmp > 1.e10, tmp)
        tmp = np.ma.masked_where(tmp <=0 , tmp)
        im2=plt.imshow(tmp, cmap=cm.YlOrRd, extent=(-180,180,-90,90))

        plt.scatter( lon, lat, s= 20, c='r')
        plt.scatter( lon_global[ix1], lat_global[iy1], s=10, c='b')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.15)

        cbar=plt.colorbar(im2, cax=cax)
        cbar.set_label( 'Average discharge (m3/s)', size=ssize)
        cbar.ax.tick_params(labelsize=ssize)

        pdf.savefig()
        plt.close()

        # =====
        fig = plt.figure(figsize=(ssize ,ssize*0.5),dpi=300)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        plt.title('GRDC: '+ str(number) + ' dis(grdc):' + str(int(grdc_mean)) + ' (cmf):'+str(int(np.nanmean(dis))))

        tmp = np.nanmean(outflw, axis=0)
        tmp = np.ma.masked_where(tmp > 1.e10, tmp)
        tmp = np.ma.masked_where(tmp <=0 , tmp)
        vmax = np.nanmax( tmp[iy1-12:iy1+12, ix1-12:ix1+12] )
        im2=plt.imshow(tmp, cmap=cm.YlOrRd, extent=(-180,180,-90,90), vmax=vmax)

        plt.scatter( lon, lat, s= 20, c='r', label='GRDC')
        plt.scatter( lon_global[ix1], lat_global[iy1], s=10, c='b', label='CMF')

        ax.set_xlim( lon - 3, lon + 3)
        ax.set_ylim( lat - 3, lat + 3)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.15)

        cbar=plt.colorbar(im2, cax=cax)
        cbar.set_label( 'Average discharge (m3/s)', size=ssize)
        cbar.ax.tick_params(labelsize=ssize)
        plt.legend(loc=0)

        pdf.savefig()
        plt.close()
        '''

        #print outflw.shape
        #print outflw[:,iy1,ix1]

        #print lon, lat, lon_global[ix1], lat_global[iy1]
    return wse, wse_max, wse_min, wse_max_loc, wse_min_loc
#======================================================================================================
def read_hweb(inputlist):
    station=inputlist[0]
    syear=inputlist[1]
    eyear=inputlist[2]
    egm08=inputlist[3]
    egm96=inputlist[4]
    hwb_data=hweb.HydroWeb_continous_WSE(station,syear,1,1,eyear,12,31,egm08,egm96)
    return hwb_data
#======================================================================================================
def mk_plot(df,val,unit,num,ax=None):
    ax=ax or plt.gca()
    ax.plot(df[val+'1'],df[val+'2'],"o",color="grey",linestyle='none',linewidth=0.0,marker="o",fillstyle="none",markersize=3)
    min_val=min([min(df[val+'1'].values),min(df[val+'2'].values)])
    max_val=max([max(df[val+'1'].values),max(df[val+'2'].values)])
    print (min_val,max_val)
    # ax.plot([min_val,min_val ],[max_val,max_val],color="k",linestyle='--',linewidth=0.5)
    ax.plot([-100.0,100.0],[-100.0,100.0],color="k",linestyle='--',linewidth=1.0)
    # ax.axline((0.0,0.0),slope=1.0, color='k', linestyle='--',linewidth=0.5)
    ax.set_xlim(min_val,max_val)
    ax.set_ylim(min_val,max_val)
    #=================
    if val=="RMSE":
        val1="RMSE"
    elif val=="bias":
        val1="Bias"
    elif val=="corr":
        val1="CC"
    #=================
    if unit=="none":
        ax.set_xlabel("$"+val1+"_{Original}$",fontsize=10)
        ax.set_ylabel("$"+val1+"_{Secondary}$",fontsize=10)
    else:
        ax.set_xlabel("$"+val1+"_{Original}$ $["+unit+"]$",fontsize=10)
        ax.set_ylabel("$"+val1+"_{Secondary}$ $["+unit+"]$",fontsize=10)
    ax.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[num-1]),ha="left",va="center",transform=ax.transAxes,fontsize=10)
    return 0
#======================================================================================================
# read altimetry list
mapname='glb_06min'
CaMa_dir='/cluster/data6/menaka/CaMa-Flood_v4'
runoff_folder='/home/yamadai/data/CaMa_v400_simulations/VIC_BC_3h_06min'
# /cluster/data6/menaka/AltiMaP/out/altimetry_glb_06min_20230407_secondary.txt
fname="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20230407_secondary.txt"
syear=2002
eyear=2019
nbdays = int( (datetime.date(eyear + 1, 1,1) - datetime.date(syear, 1, 1)). days)
nbyears = eyear - syear + 1
#======================================================================================================
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
#======================================================================================================
"""
df=pd.read_csv(fname, sep='\s+', header=0)
#
numbers=df['ID'].values
llats=df['lat'].values
llons=df['lon'].values
#======================================================================================================
# get CMF results for original location
xlist1 = df['ix'].values - 1
ylist1 = df['iy'].values - 1
wse_cmf1, wse_cmf_max, wse_cmf_min, wse_cmf_max_loc, wse_cmf_min_loc = read_wse_multi(xlist1, ylist1, syear, eyear, numbers, llats, llons)
print (wse_cmf1.shape)

# get CMF results for secondary location
xlist2 = df['ix2'].values - 1
ylist2 = df['iy2'].values - 1
wse_cmf2, wse_cmf_max, wse_cmf_min, wse_cmf_max_loc, wse_cmf_min_loc = read_wse_multi(xlist2, ylist2, syear, eyear, numbers, llats, llons)
print (wse_cmf2.shape)

inputlist=[[df["station"][i],syear,eyear,df["EGM08"][i],df["EGM96"][i]] for i in df.index.tolist()]
#=========================================
P=Pool(20)
wse_hwb=P.map(read_hweb,inputlist)
P.close()
P.join()

wse_hwb=np.array(wse_hwb)
print (wse_hwb.shape)

df['RMSE1']=np.array([rmse(s,o) for s,o in zip(wse_cmf1, wse_hwb)])
df['bias1']=np.array([bias(s,o) for s,o in zip(wse_cmf1, wse_hwb)])
df['corr1']=np.array([corr(s,o) for s,o in zip(wse_cmf1, wse_hwb)])

df['RMSE2']=np.array([rmse(s,o) for s,o in zip(wse_cmf2, wse_hwb)])
df['bias2']=np.array([bias(s,o) for s,o in zip(wse_cmf2, wse_hwb)])
df['corr2']=np.array([corr(s,o) for s,o in zip(wse_cmf2, wse_hwb)])

df.to_csv("Compare_secondary_VIC_BC.csv",sep=",",index=False)
print (df.head())
"""

df=pd.read_csv("Compare_secondary_VIC_BC.csv",sep=",",header=0)
#-------------------------------------------
hgt=11.69*(1.0/3.0)
wdt=8.27*(2.0/2.0)
fig=plt.figure(figsize=(wdt, hgt))
#plt.title(pname[point][0],fontsize=12)
G  = gridspec.GridSpec(ncols=2,nrows=1)
ax1 = fig.add_subplot(G[0,0])
mk_plot(df,"RMSE","m", 1, ax=ax1)
print (sum(df["RMSE1"].values <= df["RMSE2"].values)/float(len(df["RMSE1"].values)))

ax2 = fig.add_subplot(G[0,1])
mk_plot(df,"bias","m", 2, ax=ax2)
print (sum(df["bias1"].values <= df["bias2"].values)/float(len(df["bias1"].values)))

# ax3 = fig.add_subplot(G[0,2])
# mk_plot(df,"corr","none", 3, ax=ax3)

# ax1.plot(df['RMSE1'],df['RMSE2'],"o",color="grey",linestyle='none',linewidth=0.3,marker="o",fillstyle="none",markersize=2)
# min_val=min([min(df['RMSE1']),min(df['RMSE2'])])
# max1_val=max([max(df['RMSE1']),max(df['RMSE2'])])
# ax1.plot([min_val,min_val],[max_val,max_val],"--",color="k",linestyle='solid',linewidth=0.2)
# ax1.set_xlim(0,max([max(df['RMSE1']),max(df['RMSE2'])]))
# ax1.set_ylim(0,max([max(df['RMSE1']),max(df['RMSE2'])]))
# ax1.set_xlabel("$RMSE_1$ $[m]$",fontsize=10)
# ax1.set_ylabel("$RMSE_2$ $[m]$",fontsize=10)
# ax.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[num-1]),ha="left",va="center",transform=ax.transAxes,fontsize=10)

# ax2 = fig.add_subplot(G[0,2])
# ax2.plot(df['bias1'],df['bias2'],"o",color="grey",linestyle='none',linewidth=0.3,marker="o",fillstyle="none",markersize=2)
# min_val=min([min(df['bias1']),min(df['bias2'])])
# max_val=max([max(df['bias1']),max(df['bias2'])])
# ax2.plot([min_val,min_val],[max_val,max_val],"--",color="k",linestyle='solid',linewidth=0.2)
# ax2.set_xlim(0,max([max(df['bias1']),max(df['bias2'])]))
# ax2.set_ylim(0,max([max(df['bias1']),max(df['bias2'])]))
# ax2.set_xlabel("$Bias_1$ $[m]$",fontsize=10)
# ax2.set_ylabel("$Bias_2$ $[m]$",fontsize=10)

plt.tight_layout()
plt.savefig("./fig/fs03-Original_vs_secondary.png",dpi=800)
plt.savefig("./fig/fs03-Original_vs_secondary.pdf",dpi=800)
plt.savefig("./fig/fs03-Original_vs_secondary.jpg",dpi=800)