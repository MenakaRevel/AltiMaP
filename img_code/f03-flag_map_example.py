#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import matplotlib.patches as mpatches
import re
import os
import sys
import errno
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from matplotlib import colors
import matplotlib.cm as cm
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature

sys.path.append("./src")
# import read_hydroweb as hweb
import read_cgls as cgls
import read_hydrosat as hsat
import read_icesat as isat
import read_grrats as grt
from read_upstream import upstream
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
def get_data(station, tag):
    time=-9999.0
    data=-9999.0
    if tag=="HydroWeb":
        time,data=hweb.HydroWeb_WSE(station,syear=1992)
    elif tag=="HydroSat":
        time,data=hsat.HydroSat_WSE(station,syear=1992)
    elif tag=="ICESat":
        time,data=isat.ICESat_WSE(station,syear=1992)
    elif tag=="CGLS":
        time,data=cgls.CGLS_WSE(station,syear=1992)
    elif tag=="GRRATS":
        time,data=grt.GRRATS_WSE(station,syear=1992)
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
def flag_diff(flag):
    if flag == 10 or flag == 20:
        return 10
    if flag == 30 or flag == 32:
        return 20
    if flag == 31 or flag == 50:
        return 30
    if flag == 40:
        return 40
#==================================
def vec_par(LEVEL,ax=None,sup=2,w=0.005,width=0.5):
    '''
    Plotting the river network
    '''
    ax=ax or plt.gca()
    txt="tmp_%02d.txt"%(LEVEL)
    os.system("./bin/print_rivvec tmp1.txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    #print LEVEL, width#, lon1,lat1,lon2-lon1,lat2-lat1#x1[0],y1[0],x1[1]-x1[0],y1[1]-y1[0]
    # open tmp2.txt
    with open(txt,"r") as f:
        lines = f.readlines()

    #print LEVEL, width, lines, txt
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        ix = int((lon1 - west)*(1/gsize))
        iy = int((-lat1 + north)*(1/gsize))

        if rivermap[iy,ix] == 0:
            continue

        if lon1-lon2 > 180.0:
            # print (lon1, lon2)
            lon2=180.0
        elif lon2-lon1> 180.0:
            # print (lon1,lon2)
            lon2=-180.0
        #--------
        colorVal="w" #"grey"#
        #print (lon1,lon2,lat1,lat2,width)
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax)
#==================================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha)
#==================================
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
#=====================================
markers={"HydroWeb":"o","CGLS":"s","ICESat":"^","HydroSat":"X","GRRATS":"D"}
colors={"HydroWeb":"xkcd:reddy brown","CGLS":"xkcd:dark pink","ICESat":"xkcd:pinkish","HydroSat":"xkcd:light urple","GRRATS":"xkcd:tangerine"} 
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
syear=2003
eyear=2020
start=datetime.date(syear,1,1)
end=datetime.date(eyear,12,31)
days=(end-start).days + 1
#=============================
nums=[]
river=[]
pname=[]
lons =[]
lats =[]
xlist=[]
ylist=[]
leled=[]
egm08=[]
egm96=[]
llsat=[]
ldtom=[]
lflag=[]
kx1lt=[]
ky1lt=[]
kx2lt=[]
ky2lt=[]
#-------------------------------------------
# fname="./out/altimetry_"+mapname+"_test.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210618.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210701.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210706.txt"
# fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210909.txt"
fname="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210920.txt"
############################################################
f=open(fname,"r")
lines=f.readlines()
for line in lines[1::]:
    line    = filter(None,re.split(" ",line))
    #print line
    num     = line[0]
    station = line[1].strip()
    line2   = re.split("_",station)
    # print num, line2
    riv     = line2[1]
    stream  = line2[2]
    lon     = float(line[3])
    lat     = float(line[4])
    ix      = int(line[5])-1
    iy      = int(line[6])-1
    eled    = float(line[7])
    EGM08   = float(line[8])
    EGM96   = float(line[9])
    sat     = line[10].strip()
    dist    = float(line[11])
    flag    = int(line[12])
    kx1     = int(line[13])
    ky1     = int(line[14])
    kx2     = int(line[15])
    ky2     = int(line[16])
    #---
    # print (riv, flag)
    #---------------------------
    nums.append(num)
    river.append(riv)
    pname.append(station)
    lons.append(lon)
    lats.append(lat)
    xlist.append(ix)
    ylist.append(iy)
    leled.append(eled)
    egm08.append(EGM08)
    egm96.append(EGM96)
    llsat.append(sat)
    ldtom.append(dist)
    lflag.append(flag)
    kx1lt.append(kx1)
    ky1lt.append(ky1)
    kx2lt.append(kx2)
    ky2lt.append(ky2)
#=============================      
lnames=np.array(pname)
l_lons=np.array(lons)
l_lats=np.array(lats)
lflags=np.array(lflag)
# leledf=np.array(leledf)
# l_lons=np.array(l_lons)
# l_lats=np.array(l_lats)
N=float(len(lnames))
#=============================
flag_fig={}
fflag=[]
fsta=[]
flon=[]
flat=[]
fkx1=[]
fky1=[]
fkx2=[]
fky2=[]
# fpox=[0,0,0,0,1,2,3]
# fpoy=[0,1,2,3,3,3,3]
fpox=[0,1,2,3]
fpoy=[3,3,3,3]
fname="/cluster/data6/menaka/Altimetry/out/flag_fig.txt"
with open(fname,"r") as f:
    lines=f.readlines()
for line in lines:
    line    = filter(None,re.split(" ",line))
    flag    = line[0]
    station = line[1].strip()
    fflag.append(flag)
    fsta.append(station)
    # pname
    index = [i for i,x in enumerate(pname) if x == station][0]
    # print (station, "index: ",index)
    flon.append(lons[index])
    flat.append(lats[index])
    fkx1.append(kx1lt[index])
    fky1.append(ky1lt[index])
    fkx2.append(kx2lt[index])
    fky2.append(ky2lt[index])
    # fpox.append()
    # fpoy.append()
#=============================
mkdir("./fig")
# mkdir("./fig/criteria")
dataname="HydroWeb"
odir="/cluster/data6/menaka/Altimetry/results"
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
restag="3sec"
res=1.0/1200.0
nx =12000
ny =12000
hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
#=============================
# river width
sup=2
w=0.02
alpha=1
width=0.5

land="#C0C0C0"
water="#FFFFFF"

west=-180.0
east=180.0
north=90.0
south=-58.0

lllat = -58.
urlat = 90.
lllon = -180.
urlon = 180.

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

#cmap=make_colormap(colors_list)
#cmap=mbar.colormap("H02")
# cmap=cm.seismic
# cmap=cm.rainbow_r
# #cmap.set_under("w",alpha=0)
# cmapL=cmap #cm.get_cmap("rainbow_r")
vmin=1.0
vmax=8.0
norm=Normalize(vmin=vmin,vmax=vmax)

# bounds=np.arange(0.0,8.0,1.0)
# bounds=np.arange(0.0,6.0,1.0)
bounds=np.arange(0.0,4.0,1.0)
#-----------
# flag identity
# 10 = location was directly found
# 20 = location was on the unit-catchment outlet
# 30 = found the nearest river
# 31 = found the nearest main river
# 40 = correction for ocean grids
# 50 = bifurcation location
#cmap=colors.ListedColormap(['grey',"xkcd:ultramarine",'xkcd:clear blue','xkcd:jungle green',"xkcd:shamrock","xkcd:electric green","xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
#cmap=colors.ListedColormap(['grey','xkcd:jungle green',"xkcd:shamrock","xkcd:electric green","xkcd:ultramarine",'xkcd:clear blue',"xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
#cmap=colors.ListedColormap(['grey','xkcd:jungle green',"xkcd:shamrock","xkcd:greeny blue","xkcd:ultramarine",'xkcd:clear blue',"xkcd:sunny yellow","xkcd:neon red","xkcd:black"])
# cmap=colors.ListedColormap(['grey',"xkcd:dark seafoam",'xkcd:deep teal',"xkcd:saffron","xkcd:purpleish",'xkcd:royal',"xkcd:peacock blue","xkcd:carmine","xkcd:black"])
# cmapL=matplotlib.colors.ListedColormap(['grey',"xkcd:dark seafoam",'xkcd:deep teal',"xkcd:saffron","xkcd:purpleish",'xkcd:royal',"xkcd:peacock blue","xkcd:carmine","xkcd:black"])
# marlist={10:'o', 20:'d', 30:'+', 31:'*', 32:'^', 40:'X', 50:'1'}
# corlist={10:'green', 20:'blue',30:'purple', 31:'yellow', 32:'xkcd:lavender', 40:'red',50:'xkcd:deep teal'}
# cmapL = matplotlib.colors.ListedColormap(['green', 'blue', 'purple', 'yellow','xkcd:lavender', 'red', 'xkcd:deep teal'])

marlist={10:'o', 20:'d', 30:'s',40:'^', 50:'P'}
colors=["#4c72b0","#de8452","#55a868","#c54e51"]
corlist={10:'#4c72b0', 20:'#de8452', 30:'#55a868', 40:'#c54e51'} #,32:'xkcd:lavender',50:'xkcd:deep teal',30:'purple'}
cmapL = matplotlib.colors.ListedColormap(colors)

# cmapL.set_under("none") #"#000000",alpha=0)
# cmapL.set_over("none")
# cmapL.colorbar_extend="neither"
norml=BoundaryNorm(bounds,cmapL.N) #len(bounds)-1)
############################################################
vmin=1.0
vmax=26.0
normM=Normalize(vmin=1,vmax=26)
boundsM=np.arange(-0.5,26,1.0)
cmapM = matplotlib.colors.ListedColormap(['w','w','grey','k','w','k','w','y','w','w','blue','w','w','w','w','w','w','w','w','w','red', 'w','w','w','w','red'])
cmapM.set_under("none") #"#000000",alpha=0)
cmapM.set_over("none")
cmapM.colorbar_extend="neither"
normm=BoundaryNorm(boundsM,cmapM.N) #len(bounds)-1)
############################################################
hgt= 11.69*(1.0/2.0)
wdt= 8.27
fig= plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(4,4)
ax = fig.add_subplot(G[0:3,:],projection=ccrs.Robinson())
#-----------------------------  
ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
#
## Draw the river network ##
box="%f %f %f %f"%(west,east,north,south) 
os.system("./bin/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > tmp1.txt") 
# map(vec_par,np.arange(5,10+1,1))
map(vec_par,np.arange(2,10+1,1))
#
pnum=len(lnames)
for point in np.arange(pnum):
    # eled=leledf[point]
    name = lnames[point]
    lon  = l_lons[point]
    lat  = l_lats[point]
    flag = lflags[point]
    # c=cmapL(norml(flag))
    # print (name,lon,lat,flag)
    flag=flag_diff(flag)
    c=corlist[flag]
    m=marlist[flag]
    ax.scatter(lon,lat,s=0.1,marker=m,zorder=110,edgecolors="k", facecolors=c,linewidth=0.3,transform=ccrs.PlateCarree()) #, 
#--
imL=ax.scatter([],[],c=[],cmap=cmapL,s=0.1,vmin=vmin,vmax=vmax,norm=norml) # cmap=cmap, norm=norml
imL.set_visible(False)
#cbar=M.colorbar(im,"right",size="2%")
ax.outline_patch.set_linewidth(0.0)
# Add suppmetary figures
for point in np.arange(len(fflag)):
    flag= fflag[point]
    lon = flon[point]
    lat = flat[point]
    pox = fpox[point]
    poy = fpoy[point]
    west, south = westsouth(lat,lon)
    north = south + 10.0
    east  = west + 10.0
    cname0 = cname(lat,lon)
    # get the dimesion of the map
    dec=2
    val=0.10
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
    # print (lllat, lllon, urlat, urlon)
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
    # print (cname0)
    visual=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".visual.bin"
    # print (visual)
    visual=np.fromfile(visual,np.int8).reshape(12000,12000)
    #-----------------------------
    ax0 = fig.add_subplot(G[poy,pox])
    # ax0.text(0.0,1.1,pname[point],va="center",ha="center",transform=ax0.transAxes,fontsize=14)
    flag_ch="Flag: %d"%(int(flag)) #lflag[point])
    ax0.text(0.5,0.05,flag_ch,va="center",ha="center",transform=ax0.transAxes,fontsize=6,color="w",zorder=113)
    m = Basemap(projection='cyl',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon, lat_ts=0,resolution='c',ax=ax0)
    try:
        m.arcgisimage(service=maps[0], xpixels=1500, verbose=False)
        print ("ArcGIS map")
    except:
        # Draw some map elements on the map
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
        m.drawrivers(color='blue')
        print ("Normal map")
    #m.drawcoastlines( linewidth=0.1, color='k' )
    # m.fillcontinents(color=land,lake_color=water,zorder=99)  
    # ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
    # ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
    # #
    m.drawparallels([lllat,urlat], labels = [1,0,0,0], fontsize=4,linewidth=0,zorder=102)
    m.drawmeridians([lllon,urlon], labels = [0,0,0,1], fontsize=4,linewidth=0,zorder=102)
    data = ma.masked_less_equal(visual[npix:spix,wpix:epix],2)
    im=m.imshow(data,interpolation="nearest",origin="upper",cmap=cmapM,norm=normm,zorder=110) # interpolation="nearest",origin="upper",
    # print (lon,lat)
    # m.scatter(lon,lat,s=0.5,marker="o",zorder=110,edgecolors="g", facecolors="g")#,transform=ccrs.PlateCarree()) #, 
    ax0.plot(lon ,lat ,color="g",marker="o",markersize=2,zorder=111) #fillstyle="none",
    #================
    kx1= fkx1[point]
    ky1= fky1[point]
    lat1 = south + 10.0 - res/2.0 - ky1*res  
    lon1 = west + res/2.0 + kx1*res
    ax0.plot(lon1 ,lat1 ,color="r",marker="o",markersize=2,zorder=112) #fillstyle="none",
    #================
    kx2= fkx2[point]
    ky2= fky2[point]
    if kx2 != -9999 and ky2 != -9999:
        lat2 = south + 10.0 - res/2.0 - ky2*res  
        lon2 = west + res/2.0 + kx2*res
        ax0.plot(lon2 ,lat2 ,color="xkcd:orange",marker="*",markersize=2,zorder=112)
# colorbar
# cax=fig.add_axes([0.4,0.38,0.45,.01])
# cbar=plt.colorbar(imL,orientation="horizontal",ticks=np.arange(0.5,5.0+1.0,1.0),cax=cax) #[10,20,30,31,40,50],extend='both',ticks=np.arange(0.0,1.0+0.001,0.1) extend='both',
# cbar.set_ticklabels(['10', '20', '30', '40'])
# # cbar.set_ticklabels(['10', '20', '30', '31', '32', '40', '50'])
# # cbar.ax.set_ticklabels(['10', '20', '30', '31', '40', '50'])  # vertically oriented colorbar
# cbar.ax.tick_params(labelsize=6)
# cbar.set_label("Allocation Flags",fontsize=8)
# legend
features=[]
for i,flag in enumerate([10,20,30,40]):
    features.append(mlines.Line2D([], [], color=corlist[flag], marker=marlist[flag],
                          markersize=5, label='%d'%(flag),linewidth=0.0))
l,b,w,h=0.4,0.38,0.45,0.01
legend=plt.legend(handles=features,bbox_to_anchor=(l,b), loc="lower right",
           bbox_transform=fig.transFigure, ncol=1,  borderaxespad=0.0, frameon=False)#

plt.savefig("./fig/f03-allocation_flag_map_example.png",dpi=800)
plt.savefig("./fig/f03-allocation_flag_map_example.pdf",dpi=800)
plt.savefig("./fig/f03-allocation_flag_map_example.jpg",dpi=800)