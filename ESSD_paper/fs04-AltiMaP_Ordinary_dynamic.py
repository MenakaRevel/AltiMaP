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
#======================================================================================================
def mk_plot(df,val,hue,style,unit,num,ax=None):
    ax=ax or plt.gca()
    sns.scatterplot(x=val+'1',y=val+'2',data=df,ax=ax,hue=hue,style=style,
    palette=sns.color_palette("deep", 2),s=15,marker="o",alpha=0.5) #palette=palette,style_order=style_order,
    # ax.plot(df1[val+'1'],df2[val+'1'],"o",color="grey",linestyle='none',linewidth=0.0,marker="o",fillstyle="none",markersize=3)
    min_val=min([min(df[val+'1'].values),min(df[val+'2'].values)])
    max_val=max([max(df[val+'1'].values),max(df[val+'2'].values)])
    print (min_val,max_val)
    # ax.plot([min_val,min_val ],[max_val,max_val],color="k",linestyle='--',linewidth=0.5)
    ax.plot([min_val,max_val ],[min_val,max_val],color="k",linestyle='--',linewidth=0.5)
    # ax.plot([-100.0,100.0],[-100.0,100.0],color="k",linestyle='--',linewidth=1.0)
    # ax.axline((0.0,0.0),slope=1.0, color='k', linestyle='--',linewidth=0.5)
    min_val=max(-100.0,min_val)
    max_val=min(100.0,max_val)
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
        ax.set_xlabel("$"+val1+"_{AltiMap}$",fontsize=10)
        ax.set_ylabel("$"+val1+"_{Ordinary}$",fontsize=10)
    else:
        ax.set_xlabel("$"+val1+"_{AltiMap}$ $["+unit+"]$",fontsize=10)
        ax.set_ylabel("$"+val1+"_{Ordinary}$ $["+unit+"]$",fontsize=10)
    ax.text(-0.05,1.05,"%s)"%(string.ascii_lowercase[num-1]),ha="left",va="center",transform=ax.transAxes,fontsize=10)
    return 0
#======================================================================================================

#==============================
CaMa_dir = "/cluster/data6/menaka/CaMa-Flood_v4"
mapname= "glb_06min"
#==============================
df1=pd.read_csv('./AltiMap_comaprison_VIC_BC.csv',sep=',')
df2=pd.read_csv('./Ordinary_comaprison_VIC_BC.csv',sep=',')
df3=pd.read_csv('/cluster/data6/menaka/AltiMaP/out/biased_removed_altimetry_glb_06min_20230804.txt',sep='\s+')
df4=pd.read_csv('/cluster/data6/menaka/AltiMaP/out/biased_removed_altimetry_glb_06min_20230407.txt',sep='\s+')

df2=df2[df2['pnames'].isin(df1['pnames'])]
df3=df3[df3['station'].isin(df1['pnames'])]
df4=df4[df4['station'].isin(df1['pnames'])]
df=df3
df["RMSE1"]=df1["RMSE(AltiMaP)"]
df["bias1"]=df1["bias(AltiMaP)"]
df["RMSE2"]=df2["RMSE(ordinary)"]
df["bias2"]=df2["bias(ordinary)"]
df["Bias_flag1"]=np.array(["biased(Static)" if flag >= 900 else "Non-biased" for flag in df4["flag"].values])
df["color_order"]=np.array(["r" if flag >= 900 else "g" for flag in df4["flag"].values])
df["Bias_flag2"]=np.array(["biased(Dynamic)" if flag >= 800 else "Non-biased" for flag in df3["flag"].values])
df["style_order"]=np.array(["D" if flag >= 900 else "o" for flag in df3["flag"].values])

print (df.head())
print (df["color_order"].values)
#-------------------------------------------
hgt=11.69*(1.0/3.0)
wdt=8.27*(2.0/2.0)
fig=plt.figure(figsize=(wdt, hgt))
#plt.title(pname[point][0],fontsize=12)
G  = gridspec.GridSpec(ncols=2,nrows=1)
ax1 = fig.add_subplot(G[0,0])
mk_plot(df,"RMSE","Bias_flag1","Bias_flag2","m", 1, ax=ax1)
print (sum(df["RMSE1"].values <= df["RMSE2"].values)/float(len(df["RMSE1"].values)))

ax2 = fig.add_subplot(G[0,1])
mk_plot(df,"bias","Bias_flag1","Bias_flag2","m", 2, ax=ax2)
print (sum(df["bias1"].values <= df["bias2"].values)/float(len(df["bias1"].values)))

plt.tight_layout()
plt.savefig("./fig/fs04-AltiMaP_vs_Ordinary_biased_method.png",dpi=800)
plt.savefig("./fig/fs04-AltiMaP_vs_Ordinary_biased_method.pdf",dpi=800)
plt.savefig("./fig/fs04-AltiMaP_vs_Ordinary_biased_method.jpg",dpi=800)