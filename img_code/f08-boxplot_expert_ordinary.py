#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import sys
import os
import errno
from numpy import ma 
import matplotlib.gridspec as gridspec
import string
from scipy.fftpack import fft, ifft, fftfreq
#from slacker import Slacker
from multiprocessing import Pool
from multiprocessing import Process
import xarray as xr
from scipy.stats.kde import gaussian_kde
from matplotlib.mlab import find
import pandas as pd 
import seaborn as sns
import math
import my_colorbar as mbar
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import re
import xarray as xr
from matplotlib.mlab import find
import seaborn as sns
#
# from read_patchMS import upstream
#==============================
CaMa_dir = "/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname= "glb_06min"
#==============================
#=================================
def vec_par(LEVEL,ax=None):
    # river width
    sup=2
    w=0.02
    width=0.5
    ax=ax or plt.gca()
    txt="tmp_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec tmp1.txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    f = open(txt,"r")
    lines = f.readlines()
    f.close()
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        ix = int((lon1 + 180.)*4.0)
        iy = int((-lat1 + 90.)*4.0)

        if lon1-lon2 > 180.0:
            print lon1, lon2
            lon2=180.0
        elif lon2-lon1> 180.0:
            print lon1,lon2
            lon2=-180.0

        # if rivermap[iy,ix] == 1.0:
        colorVal="w" 
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#==============================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,transform=ccrs.PlateCarree(),linewidth=width,zorder=105,alpha=alpha)
#==============================
def mk_fig(sfcelv_rmse,pnum,lons,lats,cmap,ax=None):
    #=======================================
    # make figure
    #---------------------------------------
    print "plot map"
    ax=ax or plt.gca()
    # ax=fig.add_subplot(G[0,0],projection=ccrs.Robinson())
    ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
    #####
    #--
    box="%f %f %f %f"%(lllon,urlon,urlat,lllat) 
    # os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+"  > tmp1.txt") 
    # #map(vec_par,np.arange(1,10+1,1))
    # map(vec_par,np.arange(7,10+1,1))
    # M = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c',ax1=ax1)
    # #m.drawcoastlines( linewidth=0.1, color='k' )
    # M.fillcontinents(color=land,lake_color=water)
    # M.drawmapboundary(fill_color=water)
    # M.drawparallels(np.arange(south-2.0,north+0.1,20), labels = [1,0,0,0], fontsize=10,linewidth=0)
    # M.drawmeridians(np.arange(west,east+0.1,40), labels = [0,0,0,1], fontsize=10,linewidth=0)
    for point in np.arange(pnum):
        lon=lons[point]
        lat=lats[point]
        c=cmap(norm(sfcelv_rmse[point]))
        #print lon,lat,pname[point][0],mean_bias
        ax.scatter(lon,lat,s=0.1,marker="o",zorder=110,edgecolors=c, facecolors=c,transform=ccrs.PlateCarree())
    ax.outline_patch.set_linewidth(0.0)
    # # plt.gca().outline_patch.set_linewidth(0.0)
    # #--
    # im=ax.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm)#
    # im.set_visible(False)
    # #cbar=M.colorbar(im,"right",size="2%")
    # ax.outline_patch.set_linewidth(0.0)
    # # plt.gca().outline_patch.set_linewidth(0.0)
    # #colorbar
    # cax=fig.add_axes([0.40,0.20,0.4,.01])
    # cbar=plt.colorbar(im,orientation="horizontal",extend='max',ticks=np.arange(vmin,vmax+0.1,2.0),cax=cax) #,extend='both',ticks=np.arange(0.0,1.0+0.001,0.1)
    # cbar.ax.tick_params(labelsize=6)
    # cbar.set_label("$RMSE$ $(m)$",fontsize=8)
    return 0
#==============================
def mk_fig_div(sfcelv_rmse,pnum,lons,lats,cmap,ax=None):
    #=======================================
    # make figure
    #---------------------------------------
    print "plot map"
    ax=ax or plt.gca()
    # ax=fig.add_subplot(G[0,0],projection=ccrs.Robinson())
    ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
    #####
    #--
    box="%f %f %f %f"%(lllon,urlon,urlat,lllat) 
    # os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+"  > tmp1.txt") 
    #map(vec_par,np.arange(1,10+1,1))
    # map(vec_par,np.arange(7,10+1,1))
    # M = Basemap(projection='cyl',llcrnrlat=south,urcrnrlat=north,llcrnrlon=west,urcrnrlon=east, lat_ts=0,resolution='c',ax1=ax1)
    # #m.drawcoastlines( linewidth=0.1, color='k' )
    # M.fillcontinents(color=land,lake_color=water)
    # M.drawmapboundary(fill_color=water)
    # M.drawparallels(np.arange(south-2.0,north+0.1,20), labels = [1,0,0,0], fontsize=10,linewidth=0)
    # M.drawmeridians(np.arange(west,east+0.1,40), labels = [0,0,0,1], fontsize=10,linewidth=0)
    for point in np.arange(pnum):
        lon=lons[point]
        lat=lats[point]
        c=cmap(norm(sfcelv_rmse[point]))
        #print lon,lat,pname[point][0],mean_bias
        ax.scatter(lon,lat,s=0.1,marker="o",zorder=110,edgecolors=c, facecolors=c,transform=ccrs.PlateCarree())
    ax.outline_patch.set_linewidth(0.0)
    # # plt.gca().outline_patch.set_linewidth(0.0)
    #--
    # im=ax.scatter([],[],c=[],cmap=cmap,s=0.1,vmin=vmin,vmax=vmax,norm=norm)#
    # im.set_visible(False)
    # #cbar=M.colorbar(im,"right",size="2%")
    
    # #colorbar
    # cax=fig.add_axes([0.40,0.20,0.4,.01])
    # cbar=plt.colorbar(im,orientation="horizontal",extend='both',ticks=np.arange(vmin,vmax+0.1,2.0),cax=cax) #,extend='both',ticks=np.arange(0.0,1.0+0.001,0.1)
    # cbar.ax.tick_params(labelsize=6)
    # cbar.set_label("$\Delta$$RMSE$ $(m)$",fontsize=8)
    return 0
#====
def mk_boxplot(sfcelv_rmse1,sfcelv_rmse2,ax=None):
    ax=ax or plt.gca()
    flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
    boxprops = dict(color='grey')#facecolor='none'
    whiskerprops = dict(color='grey',linestyle="--")
    capprops = dict(color='grey')
    medianprops = dict(color='grey',linestyle="-",linewidth=1.0)
    meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
    # make boplot
    box=sns.boxplot(ax=ax,data=[sfcelv_rmse1,sfcelv_rmse2], fliersize=0.0, palette=["xkcd:coral","xkcd:teal"], whis=1.5\
        ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
        ,meanprops=meanprops,capprops=capprops,medianprops=medianprops) #"Paired"
    ax.set_xticklabels(["expert","ordinary"])
    # ax.set_xlabel(["expert","ordinary"])
    ax.set_ylim(ymin=-0.2,ymax=60.2)
    return 0
#====
# sfcelv
# odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
odir1 = '/cluster/data6/menaka/Altimetry/results/HydroWeb'
fname1 = odir1+"/hydroweb_cmf_daily_wse_VIC_BC.nc"
nc1 = xr.open_dataset(fname1)
sfcelv_hydroweb1=nc1.sfcelv_hydroweb.values
sfcelv_cmf1=nc1.sfcelv_cmf.values
lons1=nc1.lon.values
lats1=nc1.lat.values
pnames1=nc1.name.values
flags1=nc1.flag.values
pnum1=len(pnames1)
print np.shape(sfcelv_cmf1), pnum1
#- masked out -9999 values
sfcelv_cmf1=ma.masked_where(sfcelv_hydroweb1==-9999.0,sfcelv_cmf1).filled(-9999.0)
#sfcelv_hydroweb=ma.masked_where(sfcelv_hydroweb==-9999.0,sfcelv_hydroweb)
sfcelv_diff1=ma.masked_where(sfcelv_hydroweb1==-9999.0,(sfcelv_cmf1-sfcelv_hydroweb1)**2).filled(-9999.0)
sfcelv_rmse1=np.mean(ma.masked_less_equal(sfcelv_diff1,0.0),axis=0)#.compressed()#
sfcelv_rmse1=sfcelv_rmse1.filled()
print np.shape(sfcelv_rmse1), type(sfcelv_rmse1)
print sfcelv_rmse1#[0:10]
# sfcelv_bias_com=ma.masked_equal(sfcelv_bias,-9999.0).compressed()


# ====
# sfcelv
# odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
odir2 = '/cluster/data6/menaka/Altimetry/results/HydroWeb'
fname2 = odir2+"/hydroweb_cmf_daily_wse_VIC_BC_ordinary.nc"
nc2 = xr.open_dataset(fname2)
sfcelv_hydroweb2=nc2.sfcelv_hydroweb.values
sfcelv_cmf2=nc2.sfcelv_cmf.values
lons2=nc2.lon.values
lats2=nc2.lat.values
pnames2=nc2.name.values
flags2=nc2.flag.values
pnum2=len(pnames2)
print np.shape(sfcelv_cmf2), pnum2
#- masked out -9999 values
sfcelv_cmf2=ma.masked_where(sfcelv_hydroweb2==-9999.0,sfcelv_cmf2).filled(-9999.0)
#sfcelv_hydroweb=ma.masked_where(sfcelv_hydroweb==-9999.0,sfcelv_hydroweb)
sfcelv_diff2=ma.masked_where(sfcelv_hydroweb2==-9999.0,(sfcelv_cmf2-sfcelv_hydroweb2)**2).filled(-9999.0)
sfcelv_rmse2=np.mean(ma.masked_equal(sfcelv_diff2,-9999.0),axis=0)#.compressed()#
sfcelv_rmse2=sfcelv_rmse2.filled()
print np.shape(sfcelv_rmse2)
print sfcelv_rmse2[0:10]
# sfcelv_bias_com=ma.masked_equal(sfcelv_bias1,-9999.0).compressed()


os.system("mkdir -p ./fig")


hgt=11.69#*(4.0/15.0)
wdt=8.27
fig=plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(3,4)
# boxplot
ax1=fig.add_subplot(G[0,1:3])
#-----------------------------
# for i in sfcelv_rmse:
#     print i
# ax.hist(sfcelv_rmse, bins=10)
box=sns.boxplot(ax=ax1,data=[sfcelv_rmse1,sfcelv_rmse2], fliersize=0.0, palette=["xkcd:coral","xkcd:teal"], whis=1.5\
    ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
    ,meanprops=meanprops,capprops=capprops,medianprops=medianprops) #"Paired"
# sns.boxplot(sfcelv_rmse)
print np.mean(sfcelv_rmse1), np.median(sfcelv_rmse1)
print np.mean(sfcelv_rmse2), np.median(sfcelv_rmse2)
# sns.histplot(sfcelv_rmse)
# sns.distplot(ma.masked_greater(sfcelv_rmse,100.0),ax=ax, hist = True, kde = False,
#     kde_kws = {'linewidth': 1,'linestyle':'-'},bins=500,
#     label = "RMSE", color="xkcd:coral",norm_hist=True)
# ax.set_xlim(xmin=-0.2,xmax=500.2)
ax1.set_xticklabels(["expert","ordinary"])
ax1.set_xlabel(["expert","ordinary"])
ax1.set_ylim(ymin=-0.2,ymax=60.2)