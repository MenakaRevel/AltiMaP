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
import warnings; warnings.simplefilter('ignore')
#
# from read_patchMS import upstream
#==============================
CaMa_dir = "/cluster/data6/menaka/CaMa-Flood_v4"
mapname= "glb_06min"
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
def mk_boxplot_old(sfcelv_rmse1,sfcelv_rmse2,ax=None):
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
    ax.set_xticklabels(["AltiMaP","ordinary"])
    ax.set_ylabel('RMSE $(m)$', color='k',fontsize=8)
    # ax.set_xlabel(["expert","ordinary"])
    ax.set_ylim(ymin=-1.2,ymax=12.2)
    # print (ylabel)
    print ("1, mean, median", np.mean(sfcelv_rmse1), np.median(sfcelv_rmse1))
    print ("2, mean, median", np.mean(sfcelv_rmse2), np.median(sfcelv_rmse2))
    return 0
#====
def mk_boxplot(var1,var2,ylabel,num,colorlist=["xkcd:coral","xkcd:teal"],labellist=["AltiMaP","ordinary"],ylim=[0.0,2.2] ,ax=None):
    ax=ax or plt.gca()
    flierprops = dict(marker='o', markerfacecolor='none', markersize=8,linestyle='none', markeredgecolor='k')
    boxprops = dict(color='grey')#facecolor='none'
    whiskerprops = dict(color='grey',linestyle="--")
    capprops = dict(color='grey')
    medianprops = dict(color='grey',linestyle="-",linewidth=1.0)
    meanprops = dict(marker='D', markeredgecolor='black',markerfacecolor='green',markersize=8)
    # make boplot
    box=sns.boxplot(ax=ax,data=[var1,var2], fliersize=0.0, palette=colorlist, whis=1.5\
        ,meanline=True, width=0.8, linewidth=0.3, dodge=True\
        ,meanprops=meanprops,capprops=capprops,medianprops=medianprops) #"Paired"
    ax.set_xticklabels(labellist)
    ax.set_ylabel(ylabel, color='k',fontsize=8) #'RMSE $(m)$'
    # ax.set_xlabel(["expert","ordinary"])
    ax.set_ylim(ymin=ylim[0],ymax=ylim[1])
    ax.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[num-1]),ha="left",va="center",transform=ax.transAxes,fontsize=10)
    print (ylabel)
    print (labellist[0], " mean :", np.mean(var1), "median :",  np.median(var1))
    print (labellist[1], " mean :", np.mean(var2), "median :" ,np.median(var2))
    return 0
#====
def print_stat(sfcelv_rmse1,sfcelv_rmse2):
    print ("***************************************************************")
    print ("1, mean, median", np.mean(sfcelv_rmse1), np.median(sfcelv_rmse1))
    print ("2, mean, median", np.mean(sfcelv_rmse2), np.median(sfcelv_rmse2))
    return 0
# # #====
# # # sfcelv
# # # odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
# # # odir1 = '/cluster/data6/menaka/AltiMaP/results/HydroWeb_20230407'
# # # odir1 = '/cluster/data6/menaka/AltiMaP/results/HydroWeb_old'
# # odir1 = '/cluster/data6/menaka/AltiMaP/results/HydroWeb'
# # fname1 = odir1+"/hydroweb_cmf_daily_wse_VIC_BC.nc"
# # nc1 = xr.open_dataset(fname1)
# # sfcelv_hydroweb1=nc1.sfcelv_hydroweb.values
# # sfcelv_cmf1=nc1.sfcelv_cmf.values
# # lons1=nc1.lon.values
# # lats1=nc1.lat.values
# # pnames1=nc1.name.values
# # flags1=nc1.flag.values
# # pnum1=len(pnames1)
# # print np.shape(sfcelv_cmf1), pnum1
# # #- masked out -9999 values
# # sfcelv_cmf1=ma.masked_where(sfcelv_hydroweb1==-9999.0,sfcelv_cmf1).filled(-9999.0)
# # #sfcelv_hydroweb=ma.masked_where(sfcelv_hydroweb==-9999.0,sfcelv_hydroweb)
# # # sfcelv_diff1=ma.masked_where(sfcelv_hydroweb1==-9999.0,(sfcelv_cmf1-sfcelv_hydroweb1)**2).filled(-9999.0)
# # # sfcelv_rmse1=np.mean(ma.masked_less_equal(sfcelv_diff1,0.0),axis=0)#.compressed()#
# # # sfcelv_rmse1=np.sqrt(sfcelv_rmse1)
# # # sfcelv_rmse1=sfcelv_rmse1.filled()
# # # print np.shape(sfcelv_rmse1), type(sfcelv_rmse1)
# # # print sfcelv_rmse1#[0:10]
# # # sfcelv_bias_com=ma.masked_equal(sfcelv_bias,-9999.0).compressed()

# # # prepare as pandas dataframe
# # print (np.shape(pnames1), np.shape(lons1), np.shape(lats1), np.shape(flags1), np.shape(nc1.sfcelv_cmf.values))
# # print (len(pnames1), len(lons1), len(lats1), len(flags1), len(np.transpose(nc1.sfcelv_cmf.values)))
# # df1=pd.DataFrame({'pnames':pnames1.flatten(),'lons':lons1,'lats':lats1})
# # df1['flag']=[int(math.floor(flag/10.0)*10.0) for flag in flags1]
# # # for AltiMaP
# # df1['RMSE(AltiMaP)']=np.array([rmse(s,o) for s,o in zip(np.transpose(nc1.sfcelv_cmf.values), np.transpose(nc1.sfcelv_hydroweb.values))])
# # df1['bias(AltiMaP)']=np.array([bias(s,o) for s,o in zip(np.transpose(nc1.sfcelv_cmf.values), np.transpose(nc1.sfcelv_hydroweb.values))])
# # df1['corr(AltiMaP)']=np.array([corr(s,o) for s,o in zip(np.transpose(nc1.sfcelv_cmf.values), np.transpose(nc1.sfcelv_hydroweb.values))])

# # # Ordinary
# # # ====
# # # sfcelv
# # # odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
# # # odir2 = '/cluster/data6/menaka/AltiMaP/results/HydroWeb_20230407'
# # # odir2 = '/cluster/data6/menaka/AltiMaP/results/HydroWeb_old'
# # odir2 = '/cluster/data6/menaka/AltiMaP/results/HydroWeb'
# # fname2 = odir2+"/hydroweb_cmf_daily_wse_VIC_BC_ordinary.nc"
# # nc2 = xr.open_dataset(fname2)
# # sfcelv_hydroweb2=nc2.sfcelv_hydroweb.values
# # sfcelv_cmf2=nc2.sfcelv_cmf.values
# # lons2=nc2.lon.values
# # lats2=nc2.lat.values
# # pnames2=nc2.name.values
# # flags2=nc2.flag.values
# # pnum2=len(pnames2)
# # print np.shape(sfcelv_cmf2), pnum2
# # #- masked out -9999 values

# # # sfcelv_cmf2=ma.masked_where(sfcelv_hydroweb2==-9999.0,sfcelv_cmf2).filled(-9999.0)
# # # #sfcelv_hydroweb=ma.masked_where(sfcelv_hydroweb==-9999.0,sfcelv_hydroweb)
# # # sfcelv_diff2=ma.masked_where(sfcelv_hydroweb2==-9999.0,(sfcelv_cmf2-sfcelv_hydroweb2)**2).filled(-9999.0)
# # # sfcelv_rmse2=np.mean(ma.masked_equal(sfcelv_diff2,-9999.0),axis=0)#.compressed()#
# # # sfcelv_rmse2=np.sqrt(sfcelv_rmse2)
# # # sfcelv_rmse2=sfcelv_rmse2.filled()
# # # print np.shape(sfcelv_rmse2)
# # # print sfcelv_rmse2[0:10]
# # # sfcelv_bias_com=ma.masked_equal(sfcelv_bias1,-9999.0).compressed()

# # # prepare as pandas dataframe
# # df2=pd.DataFrame({'pnames':pnames2.flatten(),'lons':lons2,'lats':lats2})
# # df2['flag']=[int(math.floor(flag/10.0)*10.0) for flag in flags2]
# # # for ordinary
# # df2['RMSE(ordinary)']=np.array([rmse(s,o) for s,o in zip(np.transpose(nc2.sfcelv_cmf.values), np.transpose(nc2.sfcelv_hydroweb.values))])
# # df2['bias(ordinary)']=np.array([bias(s,o) for s,o in zip(np.transpose(nc2.sfcelv_cmf.values), np.transpose(nc2.sfcelv_hydroweb.values))])
# # df2['corr(ordinary)']=np.array([corr(s,o) for s,o in zip(np.transpose(nc2.sfcelv_cmf.values), np.transpose(nc2.sfcelv_hydroweb.values))])



# # # print (np.argwhere((sfcelv_rmse2>8.0)*flags2==10)*1 == 1)
# # # print (pnames2[0], sfcelv_rmse2[0])
# # # print (pnames2[np.argwhere((sfcelv_rmse2>10)*flags2==10)*1 == 1])
# # # # for num in np.arange(pnum2):
# # # #     if sfcelv_rmse2[num]-sfcelv_rmse1[num] > 10.0:
# # # #         print (pnames2[num][0], sfcelv_rmse2[num], sfcelv_rmse1[num])

# # # # os.system("mkdir -p ./fig")

# # print (df1.head())
# # print (df1.describe())
# # print (df1[df1['flag'] <= 200].describe())

# # print (df2.head())
# # print (df2.describe())

# # df1.to_csv('./AltiMap_comaprison_VIC_BC.csv',sep=',',index=False)
# # df2.to_csv('./Ordinary_comaprison_VIC_BC.csv',sep=',',index=False)


df1=pd.read_csv('./AltiMap_comaprison_VIC_BC.csv',sep=',')
df2=pd.read_csv('./Ordinary_comaprison_VIC_BC.csv',sep=',')

df2=df2[df2['pnames'].isin(df1['pnames'])]


# print (df1['RMSE(AltiMaP)'] - df2['RMSE(ordinary)'])

for i in df2.index:
    pname = df2['pnames'][i]
    # print i , pname
    if pname not in df1['pnames'].values:
        print 'no ', pname, ' in AltiMaP'
        continue
    # print (df1[df1["pnames"]==pname]['RMSE(AltiMaP)'].values - df2[df2["pnames"]==pname]['RMSE(ordinary)'].values)
    RMSE1 = df1[df1["pnames"]==pname]['RMSE(AltiMaP)'].values
    RMSE2 = df2[df2["pnames"]==pname]['RMSE(ordinary)'].values
    if  (RMSE1 - RMSE2) > 5.0:
        print pname, df1[df1["pnames"]==pname]['flag'].values[0], RMSE1[0], RMSE2[0] #'RMSE(AltiMaP) >  RMSE(ordinary) for ', 
    # if (RMSE1 - RMSE2) < 0.0:
    #     print 'RMSE(AltiMaP) <  RMSE(ordinary) for ', pname, df1[df1["pnames"]==pname]['flag'].values
    # else:
    #     print 'RMSE(AltiMaP) =  RMSE(ordinary) for ', pname, df1[df1["pnames"]==pname]['flag'].values



"""
hgt=11.69*(1.0/4.0)
wdt=8.27*(2.0/2.0)
fig=plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(nrows=1,ncols=3) #3 
# overall # mk_boxplot(sfcelv_rmse1,sfcelv_rmse2,ax=ax0) 
# RMSE
ax1=fig.add_subplot(G[0,0])
mk_boxplot(df1[df1['flag'] <= 200]['RMSE(AltiMaP)'],df2['RMSE(ordinary)'],'RMSE $(m)$',1,ylim=[-10.2,10.2],ax=ax1)
print ("pandas dataframe",df1[df1['flag'] <= 200]['RMSE(AltiMaP)'].median(),df2['RMSE(ordinary)'].median())
print ("numpy array",np.median(sfcelv_rmse1), np.median(sfcelv_rmse2))
# bias
ax2=fig.add_subplot(G[0,1])
mk_boxplot(df1[df1['flag'] <= 200]['bias(AltiMaP)'],df2['bias(ordinary)'],'bias $(m)$',2, ylim=[-10.2,10.2],ax=ax2)

# RMSE
ax3=fig.add_subplot(G[0,2])
mk_boxplot(df1[df1['flag'] <= 200]['corr(AltiMaP)'],df2['corr(ordinary)'],'CC',3,ylim=[-1.2,1.2],ax=ax3)

# group by flag
print ("group by flag")
print ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print ("AltiMaP")
print ("All")
print (df1[['RMSE(AltiMaP)','bias(AltiMaP)','corr(AltiMaP)']].describe())
print ("Flag-wise")
print (df1.groupby('flag')['RMSE(AltiMaP)','bias(AltiMaP)','corr(AltiMaP)'].describe())
print ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print ("Ordinary")
print ("All")
print (df2[['RMSE(ordinary)','bias(ordinary)','corr(ordinary)']].describe())
print ("Flag-wise")
print ("Ordinary",df2.groupby('flag')['RMSE(ordinary)','bias(ordinary)','corr(ordinary)'].describe())
#===============================================================================
# print ("Flag 10")
# print_stat(sfcelv_rmse1[np.logical_or(flags1==10,flags1==11,flags1==12)],sfcelv_rmse2[np.logical_or(flags2==10,flags2==11,flags2==12)])
# print ("Flag 20")
# print_stat(sfcelv_rmse1[np.logical_or(flags1==20,flags1==21)],sfcelv_rmse2[np.logical_or(flags2==20,flags2==21)])
# print ("Flag 30")
# print_stat(sfcelv_rmse1[np.logical_or(flags1==30,flags1==31)],sfcelv_rmse2[np.logical_or(flags2==30,flags2==31)])
# print ("Flag 40")
# print_stat(sfcelv_rmse1[np.where(flags1==40)],sfcelv_rmse2[np.where(flags2==40)])
# ax0.text(-0.05,1.05,"%s) All"%(string.ascii_lowercase[0]),ha="left",va="center",transform=ax0.transAxes,fontsize=10)
# print np.shape(sfcelv_rmse1), np.shape(flags1)
# print sfcelv_rmse1[np.logical_or(flags1==10,flags1==20)]
# # Flag 10
# ax1=fig.add_subplot(G[1,0:2])
# mk_boxplot(sfcelv_rmse1[np.logical_or(flags1==10,flags1==20)],sfcelv_rmse2[np.logical_or(flags2==10,flags2==20)],ax=ax1)
# ax1.text(-0.05,1.05,"%s) Flag 10"%(string.ascii_lowercase[1]),ha="left",va="center",transform=ax1.transAxes,fontsize=10)
# # Flag 20
# ax2=fig.add_subplot(G[1,2::])
# mk_boxplot(sfcelv_rmse1[np.logical_or(flags1==30,flags1==32)],sfcelv_rmse2[np.logical_or(flags2==30,flags2==32)],ax=ax2)
# ax2.text(-0.05,1.05,"%s) Flag 20"%(string.ascii_lowercase[2]),ha="left",va="center",transform=ax2.transAxes,fontsize=10)
# # Flag 30
# ax3=fig.add_subplot(G[2,0:2])
# mk_boxplot(sfcelv_rmse1[np.logical_or(flags1==31,flags1==50)],sfcelv_rmse2[np.logical_or(flags2==31,flags2==50)],ax=ax3)
# ax3.text(-0.05,1.05,"%s) Flag 30"%(string.ascii_lowercase[3]),ha="left",va="center",transform=ax3.transAxes,fontsize=10)
# # Flag 40
# ax4=fig.add_subplot(G[2,2::])
# mk_boxplot(sfcelv_rmse1[np.where(flags1==40)],sfcelv_rmse2[np.where(flags2==40)],ax=ax4)
# ax4.text(-0.05,1.05,"%s) Flag 40"%(string.ascii_lowercase[4]),ha="left",va="center",transform=ax4.transAxes,fontsize=10)
# plt.show()
plt.tight_layout()
plt.savefig("./fig/f09-boxplot_expert_ordinary_flag.png",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./fig/f09-boxplot_expert_ordinary_flag.jpg",dpi=800,bbox_inches="tight", pad_inches=0.0)
plt.savefig("./fig/f09-boxplot_expert_ordinary_flag.pdf",dpi=800,bbox_inches="tight", pad_inches=0.0)
"""