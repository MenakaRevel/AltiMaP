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
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import re
import sys
import os
import errno
from scipy import stats
from sklearn.metrics import r2_score
import pandas as pd

sys.path.append("../src")
# from read_patchMS import upstream
# from river_function import river_profile
# import read_hydroweb as hweb
# import read_cgls as cgls
# import read_hydrosat as hsat
# import read_icesat as isat
# import read_grrats as grt
#=========================================
def get_data(station, tag, syear=2000, eyear=2020,egm08=0.0,egm96=0.0):
    time=-9999.0
    data=-9999.0
    if tag=="HydroWeb":
        time,data=hweb.HydroWeb_WSE(station,syear=syear,eyear=eyear,egm08=egm08,egm96=egm96)
    elif tag=="HydroSat":
        time,data=hsat.HydroSat_WSE(station,syear=syear,eyear=eyear,egm08=egm08,egm96=egm96)
    elif tag=="ICESat":
        time,data=isat.ICESat_WSE(station,syear=syear,eyear=eyear,egm08=egm08,egm96=egm96)
    elif tag=="CGLS":
        time,data=cgls.CGLS_WSE(station,syear=syear,eyear=eyear,egm08=egm08,egm96=egm96)
    elif tag=="GRRATS":
        time,data=grt.GRRATS_WSE(station,syear=syear,eyear=eyear,egm08=egm08,egm96=egm96)
    return time, data
#=============================
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
def river_along_profile(ix,iy,west,south,csize,nx,ny,hiresmap):
    length,elevation,k = river_profile(ix,iy,west,south,csize,nx,ny,hiresmap)
    # print (length, elevation, k)
    return length[0:k] , elevation[0:k]
#=============================
def r_squared(x,y,a,b):
    # for linear function of y = a + bx 
    # fit values, and mean
    x    = np.array(x)
    y    = np.array(y)
    yhat = a + b * x                 # or [a + b * z for z in x]
    # print (y)
    # print (yhat)
    # for i in range(len(x)):
    #     print (x[i], y[i], yhat[i])
    ybar  = np.sum(y)/len(y)            # or sum(y)/len(y) # np.mean(y) 
    sstot = np.sum((y    - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    ssreg = np.sum((yhat - ybar)**2)    # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sserr = np.sum((yhat - y   )**2)    # or sum([ (yhat[i] - y[i])**2 for i in range(len(y))])
    # print ("$r^2$ calculation: ",ssreg, sserr, sstot, ssreg+sserr)
    r2 = ssreg / sstot

    return r2
#=============================
def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
#=============================
def covert_lonlat(lon, lat, west=-180.0, north=90.0, gsize=0.1):
    '''
    Convert lat lon to x, y coordinate
    '''
    ix = int((lon - west)*(1/gsize))
    iy = int((-lat + north)*(1/gsize))

    return ix, iy
#=============================
mkdir("../fig")
mkdir("../fig/SWOT_high_res_map")
#=============================
# sfcelv
syear=2000
eyear=2020
start=datetime.date(syear,1,1)
end=datetime.date(eyear,12,31)
days=(end-start).days + 1
#=========================================
#TAG="CGLS"
TAG="HydroWeb"
# TAG="ICESat"
# TAG="HydroSat"
#========================================= 
station0=82215000310101 #82242900080501 #82219000220331 #82219000060011 #82215000220011 # 82213000040021 #82219000280291
dataname="SWOT"
odir="/cluster/data6/menaka/AltiMaP/results"
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
restag="3sec"
#=========================================
obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20241129.txt"
df=pd.read_csv(obstxt, sep='\s+', header=0) #,encoding=str)
print (df.columns)
df["flag"]=np.array([int(math.floor(flag/10.0)*10.0) for flag in df["flag"]])
# df["uparea"]=[uparea[iy-1,ix-1]*1e-6 for ix, iy in zip(df["ix"],df["iy"])]
N=float(len(df))
print (df.head())
#=============================
TAG=dataname
res=1.0/1200.0
nx =12000
ny =12000
hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
if restag == "3sec":
    res=1.0/1200.0
    nx =12000
    ny =12000
    hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
############################################################
# pnum=10 #len(pname)
#print np.shape(sfcelv_hydroweb)
colors=['xkcd:pastel blue','xkcd:teal','xkcd:aqua green','xkcd:dark pink','xkcd:purple','xkcd:magenta']
labels=["cmf oroginal","cmf interpolated","cmf ele diff",TAG]
#=============================
vmin=1.0
vmax=26.0
norm=Normalize(vmin=vmin,vmax=vmax)
bounds=np.arange(-0.5,26,1.0)
############################################################
#=====================================
# high-resolution data
# visualization: 
# visual
# 0  - sea
# 1  - land(undefied)
# 2  - land(defined in CaMa)
# 3  - grid box
# 5  - catchment boundry
# 10 - channel
# 20 - outlet pixel
# 25 - river mouth
#==============
# cmapL = matplotlib.colors.ListedColormap(['w','w','grey','k','w','k','w','y','w','w','blue','w','w','w','w','w','w','w','w','w','red', 'w','w','w','w','red'])
cmapL = matplotlib.colors.ListedColormap(['w','#044830','#044830','w','#0070FF','k','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','#0070FF','red', '#0070FF','#0070FF','#0070FF','#0070FF','red'])
cmapL.set_under("none") #"#000000",alpha=0)
cmapL.set_over("none")
cmapL.colorbar_extend="neither"
norml=BoundaryNorm(bounds,cmapL.N) #len(bounds)-1)
cmapR = matplotlib.colors.ListedColormap(['b'])
############################################################
# pnum=len(pname)
#=====================================
markers={"HydroWeb":"o","CGLS":"s","ICESat":"^","HydroSat":"X","GRRATS":"D"}
colors={"HydroWeb":"xkcd:sea blue","CGLS":"xkcd:dark pink","ICESat":"xkcd:pinkish","HydroSat":"xkcd:light urple","GRRATS":"xkcd:tangerine"} 
# reddy brown
#=============================
maps = ['ESRI_Imagery_World_2D',    # 0
        'ESRI_StreetMap_World_2D',  # 1
        'NatGeo_World_Map',         # 2
        'NGS_Topo_US_2D',           # 3
        'Ocean_Basemap',            # 4
        'USA_Topo_Maps',            # 5
        'World_Imagery',            # 6
        'World_Physical_Map',       # 7
        'World_Shaded_Relief',      # 8
        'World_Street_Map',         # 9
        'World_Terrain_Base',       # 10
        'World_Topo_Map'            # 11
        ]
#=============================
#============================
# with PdfPages(pdfname) as pdf:
# for point in np.arange(0,pnum):
print (station0)
print (df["station"].values[0:10])
pname=df[df["station"]==station0]["station"].values
lflag=df[df["station"]==station0]["flag"].values
lons=df[df["station"]==station0]["lon"].values
lats=df[df["station"]==station0]["lat"].values
kx1lt=df[df["station"]==station0]["kx1"].values
ky1lt=df[df["station"]==station0]["ky1"].values
kx2lt=df[df["station"]==station0]["kx2"].values
ky2lt=df[df["station"]==station0]["ky2"].values
print (pname)
point=0
######################
print ("=======================================")
print (point, pname[point], lflag[point])
hgt=11.69*(1.0/2.0)
wdt=8.27
fig=plt.figure(figsize=(wdt, hgt))
#plt.title(pname[point][0],fontsize=12)
# G = gridspec.GridSpec(3,2)
G = gridspec.GridSpec(1,1)
lon = lons[point]
lat = lats[point]
west, south = westsouth(lat,lon)
north = south + 10.0
east  = west + 10.0
cname0 = cname(lat,lon)
# get the dimesion of the map
dec=2
val=0.07
lllat = round_half_down(lat-val,dec)
urlat = round_half_up(lat+val,dec)
lllon = round_half_down(lon-val,dec)
urlon = round_half_up(lon+val,dec)
if abs(lllat-urlat) < val:
    urlat=round_half_up(urlat+val,dec)
    lllat=round_half_down(lllat-val,dec)
if abs(lllon-urlon) < val:
    urlon=round_half_up(urlon+val,dec)
    lllon=round_half_down(lllon-val,dec)
#---------------------
lllat=max(lllat,south)
urlat=min(urlat,north)
lllon=max(lllon,west)
urlon=min(urlon,east)
print (lllat, lllon, urlat, urlon)
#=====================================
# londiff=int((east-west)*1200)
# latdiff=int((north-south)*1200)
npix= int((north-urlat)*1200)
spix= int((north-lllat)*1200)
wpix= int((lllon-west)*1200)
epix= int((urlon-west)*1200)
print (npix,":",spix,",",wpix,":",epix)
#=====================================
# high-resolution data
# visualization: 
# visual
# 0  - sea
# 1  - land(undefied)
# 2  - land(defined in CaMa)
# 3  - grid box
# 5  - catchment boundry
# 10 - channel
# 20 - outlet pixel
# 25 - river mouth
#==============
# print (cname0)
# visual
visual=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".visual.bin"
# print (visual)
visual=np.fromfile(visual,np.int8).reshape(12000,12000)
# rivwidth
rivwth=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".rivwth.bin"
# print (rivwth)
rivwth=np.fromfile(rivwth,np.float32).reshape(12000,12000)
#-----------------------------
ax0 = fig.add_subplot(G[0,0])
# ax0.text(0.0,1.1,pname[point],va="center",ha="center",transform=ax0.transAxes,fontsize=14)
# flag_ch="flag: %d"%(lflag[point])
# ax0.text(0.8,1.1,flag_ch,va="center",ha="center",transform=ax0.transAxes,fontsize=14)
m = Basemap(projection='cyl',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon, lat_ts=0,resolution='c',ax=ax0)
try:
    # m.arcgisimage(service=maps[1], xpixels=1500, verbose=False)
    m.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=1200)
    print ("ArcGIS map")
except:
    # Draw some map elements on the map
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    # m.drawrivers(color='blue')
    print ("Normal map")
#m.drawcoastlines( linewidth=0.1, color='k' )
# m.fillcontinents(color=land,lake_color=water,zorder=99)  
# ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
# ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
# #
# m.drawparallels([lllat,urlat], labels = [1,0,0,0], fontsize=10,linewidth=0.1,zorder=102)
# m.drawmeridians([lllon,urlon], labels = [0,0,0,1], fontsize=10,linewidth=0.1,zorder=102)
data0 = ma.masked_less_equal(visual[npix:spix,wpix:epix],-9999)
im0=m.imshow(data0,interpolation="nearest",origin="upper",cmap=cmapL,norm=norml,zorder=110) # interpolation="nearest",origin="upper",
# data1 = ma.masked_greater_equal(rivwth[npix:spix,wpix:epix],0.0)
# im1=m.imshow(data1,interpolation="nearest",origin="upper",cmap=cmapR,norm=norml,zorder=110) # interpolation="nearest",origin="upper",
# print (lon,lat)
# m.scatter(lon,lat,s=0.5,marker="o",zorder=110,edgecolors="g", facecolors="g")#,transform=ccrs.PlateCarree()) #, 
# ax0.plot(lon ,lat ,color="g",marker="o",label="intitial",markersize=7,linewidth=0,zorder=111) #fillstyle="none",
#================
for point2 in range(len(df)):
    kx1= df['kx1'].values[point2]
    ky1= df['ky1'].values[point2]
    lat1 = south + 10.0 - res/2.0 - ky1*res  
    lon1 = west + res/2.0 + kx1*res
    ax0.plot(lon1 ,lat1 ,color="r",marker="o",markeredgecolor='k',label="SWOT",markersize=6,linewidth=0,zorder=112) #fillstyle="none",
#========================================================
# plt.legend(loc="upper center", bbox_to_anchor=(0.5,0.0), ncol=4)
#========================================================
plt.tight_layout()
plt.savefig("../fig/SWOT_high_res_map/SWOT_"+str(station0)+".png",dpi=500)


# # #================
# # kx2= kx2lt[point]
# # ky2= ky2lt[point]
# # if kx2 != -9999 and ky2 != -9999:
# #     lat2 = south + 10.0 - res/2.0 - ky2*res  
# #     lon2 = west + res/2.0 + kx2*res
# #     ax0.plot(lon2 ,lat2 ,color="xkcd:orange",marker="*",label="secondary",markersize=7,linewidth=0,zorder=112)
# # #================
# # # if ordinary allocation used
# # gsize=0.1 # glb_06min
# # west0=-180.0
# # north0=90.0
# # # kx3=int((lon1 - west0)/gsize)
# # # ky3=int((north0 - lat1)/gsize)
# # kx3, ky3 = covert_lonlat(lon1, lat1)
# # lat3=north0 - ky3*gsize
# # lon3=west0 + kx3*gsize
# # ax0.plot(lon3 ,lat3 ,color="xkcd:hot pink",marker="D",label="ordinary",markersize=7,linewidth=0,zorder=112)
# # # print (kx,ky,lon0,lat0)
# # # # #========================================================