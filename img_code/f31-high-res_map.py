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

sys.path.append("./src")
# from read_patchMS import upstream
from river_function import river_profile
import read_hydroweb as hweb
import read_cgls as cgls
import read_hydrosat as hsat
import read_icesat as isat
import read_grrats as grt
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
mkdir("./fig")
mkdir("./fig/high_res_map")
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
############################################################
# # rivername0=sys.argv[1] #"CONGO" #"AMAZONAS"
# # stream0=sys.argv[2]
# # dataname=sys.argv[3]
# # odir=sys.argv[4] #"/cluster/data6/menaka/AltiMaP/results"
# # mapname=sys.argv[5] #"glb_06min"
# # CaMa_dir=sys.argv[6] #"/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# # restag=sys.argv[7] #"3sec"
# # obstxt=sys.argv[8] #"./out/altimetry_"+mapname+"_test.txt"
# # stream0=["AMAZONAS","CONGO"]
#=
rivername0="CONGO" #"AMAZONAS"
stream0="CONGO" #
# station0="R_CONGO_CONGO_KM2054"
# station0="R_CONGO_CONGO_KM2053"
# station0="R_CONGO_CONGO_KM0032"
# station0="R_CONGO_CONGO_KM0031"
# station0="R_CONGO_CONGO_KM0030"
# station0="R_CONGO_BUSIRA_KM1444"
# station0="R_CONGO_CONGO_KM0309"
# station0="R_AMAZONAS_AMAZONAS_KM3491"
# station0="R_LENA_LENA_KM2195"
# station0="R_LENA_LENA_KM2501"
# station0="R_AMAZONAS_MADEIRA_KM2310"
# station0="R_AMAZONAS_MADEIRA_KM2053"
# station0="R_AMAZONAS_AMAZONAS_KM0343"
# station0="R_PARAIBA-DO-SOL_URURAI-TRIB-01_KM0109"
# station0="R_CONGO_LOMAMI_KM2608"
# station0="R_CUVO-QUEVE_CUVO-QUEVE_KM0001"
# station0="R_AMAZONAS_APAPORIS_KM2843"
# station0="R_AMAZONAS_BIA_KM2611"
# station0="R_AMAZONAS_AMAZONAS_KM0421"
# station0="R_AMAZONAS_AMAZONAS_KM1685"
# station0="R_AMAZONAS_AMAZONAS_KM0755"
# station0="R_AMAZONAS_SOLIMOES_KM1339"
# station0="R_AMAZONAS_SOLIMOES_KM1340"
# station0="R_AMAZONAS_SOLIMOES_KM1553"
# station0="R_AMAZONAS_SOLIMOES_KM1692"
# station0="R_AMAZONAS_SOLIMOES_KM1976"
# station0="R_AMAZONAS_SOLIMOES_KM2033"
# station0="R_GARONNE_GARONNE_KM0192"
# station0="R_VOLGA_UFA_KM3232"
# station0="R_VOLGA_OKA_KM3076"
station0="R_AMAZONAS_CONONACO_KM4357"
dataname="HydroWeb"
odir="/cluster/data6/menaka/AltiMaP/results"
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
restag="3sec"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210709.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210714.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210715.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210806.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210807.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210808.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210812.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210817.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210909.txt"
obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"+mapname+"_20210920.txt"
# stream0=["AMAZONAS","CONGO"]
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
# cmapL = matplotlib.colors.ListedColormap(['w','w','grey','k','w','k','w','y','w','w','blue','w','w','w','w','w','w','w','w','w','red', 'w','w','w','w','red'])
cmapL = matplotlib.colors.ListedColormap(['w','grey','grey','k','y','k','y','y','y','y','b','b','y','y','y','y','y','y','y','y','red', 'y','y','y','y','red'])
cmapL.set_under("none") #"#000000",alpha=0)
cmapL.set_over("none")
cmapL.colorbar_extend="neither"
norml=BoundaryNorm(bounds,cmapL.N) #len(bounds)-1)
cmapR = matplotlib.colors.ListedColormap(['b'])
############################################################
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
# fname="./out/altimetry_"+mapname+"_20210518.txt"
fname=obstxt
#--
with open(fname,"r") as f:
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
    #-----------------------
    # # # print (riv,station,kx,ky)
    # # if riv != rivername0:
    # #     continue
    # # if riv == "AMAZONAS":
    # #     if stream != stream0:
    # #         continue
    # # elif riv == "CONGO":
    # #     if stream != stream0:
    # #         continue
    # # else:
    # #     if stream != stream0:
    # #         continue
    # print station
    if station==station0:
        print "station found: ",station0
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
    # print (riv,station)
#=============================
pnum=len(pname)
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
# #=============================
# if TAG=="HydroWeb":
#     pdfname=odir+"/HydroWeb/high_res/hydroweb_cmf_hres_map_"+rivername0+"-"+stream0+".pdf"
# if TAG=="CGLS":
#     pdfname=odir+"/CGLS/high_res/cgls_cmf_hres_map_"+rivername0+"-"+stream0+".pdf"
# if TAG=="ICESat":
#     pdfname=odir+"/ICESat/high_res/icesat_cmf_hres_map_"+rivername0+"-"+stream0+".pdf"
# if TAG=="HydroSat":
#     pdfname=odir+"/HydroSat/high_res/hydrosat_cmf_hres_map_"+rivername0+"-"+stream0+".pdf"
# #============================
# with PdfPages(pdfname) as pdf:
# for point in np.arange(0,pnum):
print pname
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
ax0.text(0.0,1.1,pname[point],va="center",ha="center",transform=ax0.transAxes,fontsize=14)
flag_ch="flag: %d"%(lflag[point])
ax0.text(0.8,1.1,flag_ch,va="center",ha="center",transform=ax0.transAxes,fontsize=14)
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
m.drawparallels([lllat,urlat], labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
m.drawmeridians([lllon,urlon], labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
data0 = ma.masked_less_equal(visual[npix:spix,wpix:epix],-9999)
im0=m.imshow(data0,interpolation="nearest",origin="upper",cmap=cmapL,norm=norml,zorder=110) # interpolation="nearest",origin="upper",
data1 = ma.masked_greater_equal(rivwth[npix:spix,wpix:epix],0.0)
im1=m.imshow(data1,interpolation="nearest",origin="upper",cmap=cmapR,norm=norml,zorder=110) # interpolation="nearest",origin="upper",
# print (lon,lat)
# m.scatter(lon,lat,s=0.5,marker="o",zorder=110,edgecolors="g", facecolors="g")#,transform=ccrs.PlateCarree()) #, 
ax0.plot(lon ,lat ,color="g",marker="o",label="intitial",markersize=7,linewidth=0,zorder=111) #fillstyle="none",
#================
kx1= kx1lt[point]
ky1= ky1lt[point]
lat1 = south + 10.0 - res/2.0 - ky1*res  
lon1 = west + res/2.0 + kx1*res
ax0.plot(lon1 ,lat1 ,color="r",marker="o",label="best-loc",markersize=7,linewidth=0,zorder=112) #fillstyle="none",
#================
kx2= kx2lt[point]
ky2= ky2lt[point]
if kx2 != -9999 and ky2 != -9999:
    lat2 = south + 10.0 - res/2.0 - ky2*res  
    lon2 = west + res/2.0 + kx2*res
    ax0.plot(lon2 ,lat2 ,color="xkcd:orange",marker="*",label="secondary",markersize=7,linewidth=0,zorder=112)
#================
# if ordinary allocation used
gsize=0.1 # glb_06min
west0=-180.0
north0=90.0
# kx3=int((lon1 - west0)/gsize)
# ky3=int((north0 - lat1)/gsize)
kx3, ky3 = covert_lonlat(lon1, lat1)
lat3=north0 - ky3*gsize
lon3=west0 + kx3*gsize
ax0.plot(lon3 ,lat3 ,color="xkcd:hot pink",marker="D",label="ordinary",markersize=7,linewidth=0,zorder=112)
# print (kx,ky,lon0,lat0)
# # #========================================================
# # ax1 = fig.add_subplot(G[0,1])
# # length=[1.0]
# # elevation=[1.0]
# # # print ("calculate river profile...........",kx,ky,visual[ky-1,kx-1])
# # try:
# #     length, elevation = river_along_profile(kx1,ky1,west,south,res,nx,ny,hiresmap)
# # except:
# #     length=[0.0]
# #     elevation=[0.0]
# # # # print (length, elevation)
# # if visual[ky1-1,kx1-1]!=10: # or visual[ky-1,kx-1]!=20:
# #     length=[0.0]
# #     elevation=[0.0]
# # ax1.plot(length,elevation,color="grey",linestyle='-',linewidth=3.0)
# # # print (length, elevation)
# # elevtn=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".elevtn.bin"
# # # print (elevtn)
# # elevtn=np.fromfile(elevtn,np.float32).reshape(12000,12000)
# # disttom=ldtom[point]
# # # print (length[-1],"-",disttom,elevtn[ky-1,kx-1])
# # ax1.plot([length[-1]-disttom*1e3],[elevtn[ky1-1,kx1-1]],color="r",marker="o",markersize=7,linestyle='none',linewidth=0.0)
# # ax1.plot([0,length[-1]],[elevation[0],elevation[-1]],color="g",linestyle='--',linewidth=0.5)
# # locs,org = get_data(pname[point],TAG,egm08=egm08[point],egm96=egm96[point])
# # ax1.axhline(y=np.mean(org),xmin=0.0,xmax=length[-1],color=colors[TAG],linestyle='--',linewidth=0.5)
# # ax1.set_xlim(xmin=0,xmax=length[-1])
# # ax1.set_xlabel("Distance $(m)$")
# # ax1.set_ylabel("Elevation $(m EGM96)$")
# # #========================================================
# # ax2 = fig.add_subplot(G[1,:])
# # ax2.plot(locs,org,color=colors[TAG],label=TAG,linestyle='-',linewidth=0.5,marker=markers[TAG],fillstyle="none",markersize=5)
# # ax2.set_xlim(xmin=0,xmax=days+1)
# # xxlab=np.arange(syear,eyear+1,5)
# # dt=int(math.ceil(((eyear-syear)+2)/5.0))
# # xxlist=np.linspace(0,days,dt,endpoint=True)
# # ax2.set_xticks(xxlist)
# # ax2.set_xticklabels(xxlab,fontsize=8)
# # ax2.set_xlabel("Year")
# # ax2.set_ylabel("WSE $(m EGM96)$")
#========================================================
plt.legend(loc="upper center", bbox_to_anchor=(0.5,0.0), ncol=4)
#========================================================
# plt.show()
plt.savefig("./fig/high_res_map/"+station0+".png",dpi=500)
# pdf.savefig()  # saves the current figure into a pdf page
# plt.close()
# # set the file's metadata via the PdfPages object:
# d = pdf.infodict()
# d['Title'] = TAG+' high resolution map '+rivername0
# d['Author'] = 'Menaka Revel'
# d['Subject'] = TAG+' high resolution map '+rivername0
# d['Keywords'] = TAG+', high resolution, '+restag+', '+rivername0
# d['CreationDate'] = datetime.datetime(2021, 5, 25)
# d['ModDate'] = datetime.datetime.today()