#!/opt/local/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import sys
import os
import math
from numpy import ma
import itertools
from multiprocessing import Pool
from multiprocessing import Process
import re
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
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
def sat_orbit(lllon,urlon,lllat,urlat):
    with open('./Visu_TP_Tracks_HiRes_GE_V2.kml') as f:
    # with open('./cfosat_trace_v2.1.kmz') as f:
        data=f.readlines()

    incoord=0 #
    orbit={} #[track,lat,lon,index]
    lats ={}
    lons ={}
    index=0
    for line in data:
        if line.find("Ground Track")>=0:
            track=int(line[19:22])

        if line.find("</coordinates>")>=0:
            incoord=0

        if incoord==1:
            lon=float(line[:line.find(",")])
            left=line[line.find(",")+1:]
            lat=float(left[:left.find(",")])
            if lon > lllon - 25.0 and lon < urlon + 25.0 and lat > lllat - 25.0 and lat < urlat + 25.0:   
                #======================
                if track not in lons.keys():
                    lons[track]=[lon]
                else:
                    lons[track].append(lon)
                #======================
                if track not in lats.keys():
                    lats[track]=[lat]
                else:
                    lats[track].append(lat)
                #======================
                if track not in orbit.keys():
                    orbit[track]=[lat,lon,index]
                else:
                    orbit[track].append([lat,lon,index])
                index=index+1
            else:
                continue

        if line.find("<coordinates>")>=0:
            incoord=1

        if line.find("#DOT")>=0:
            break

    return orbit, lons, lats
#=============================
def vs_locations(vslist):
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
    fname="/cluster/data6/menaka/Altimetry/out/altimetry_glb_06min_20210920.txt"
    ############################################################
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        #print line
        num     = line[0]
        station = line[1].strip()
        if station not in vslist:
            continue
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
    kx1lt=np.array(kx1lt)
    ky1lt=np.array(ky1lt)
    # leledf=np.array(leledf)
    # l_lons=np.array(l_lons)
    # l_lats=np.array(l_lats)
    return lnames, l_lons, l_lats
#=============================
def vs_locations_one(vs):
    #-------------------------------------------
    fname="/cluster/data6/menaka/Altimetry/out/altimetry_glb_06min_20210920.txt"
    ############################################################
    with open(fname,"r") as f:
        lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        #print line
        num     = line[0]
        station = line[1].strip()
        if station == vs:
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
            break
    return lon, lat, kx1, ky1
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

lllat = -5.#-10.
urlat = 0.#5.
lllon = -70.#10.
urlon = -65.#25.

londiff=(east-west)*4
latdiff=(north-south)*4

npix=(90-north)*4
spix=(90-south)*4
wpix=(180+west)*4
epix=(180+east)*4

#
vslist=["R_AMAZONAS_SOLIMOES_KM2637"\
,"R_AMAZONAS_ICA_KM2615" \
,"R_AMAZONAS_JAPURA_KM2516" \
,"R_AMAZONAS_JAPURA_KM2484"\
,"R_AMAZONAS_SOLIMOES_KM2113"\
,"R_AMAZONAS_JAPURA_KM2122"\
,"R_AMAZONAS_NEGRO_KM2168"\
]

#figure
############################################################
hgt= 11.69*(1.0/3.0)
wdt= 8.27
fig= plt.figure(figsize=(wdt, hgt))
G  = gridspec.GridSpec(1,2)
ax1= fig.add_subplot(G[0,0])#,projection=ccrs.Robinson())
#-----------------------------
m = Basemap(projection='cyl',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon, lat_ts=0,resolution='c',ax=ax1)
try:
    # m.arcgisimage(service=maps[0], xpixels=1500, verbose=True)
    m.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=1200)
    print ("ArcGIS map")
except:
    # Draw some map elements on the map
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.drawrivers(color='blue')
    print ("Normal map")
# ax.set_extent([lllon,urlon,lllat,urlat],crs=ccrs.PlateCarree())
# ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor='face', facecolor=land),zorder=100)
#
orbit, lons, lats=sat_orbit(lllon,urlon,lllat,urlat)
for key in lons.keys(): #[0:10]:
    m.plot(lons[key],lats[key],color="brown",linewidth=5,alpha=0.5,zorder=110)
    # m.scatter(lons[key],lats[key])
lnames, l_lons, l_lats = vs_locations(vslist)
for lon, lat in zip(l_lons, l_lats):
    m.scatter(lon,lat,color="y",marker="o",zorder=112)
#=====================================================
ax2=fig.add_subplot(G[0,1])
############################################################
vmin=1.0
vmax=26.0
normM=Normalize(vmin=1,vmax=26)
boundsM=np.arange(-0.5,26,1.0)
cmapM = matplotlib.colors.ListedColormap(['w','w','grey','k','w','k','w','blue','w','w','blue','w','w','w','w','w','w','w','w','w','red', 'w','w','w','w','red'])
cmapM.set_under("none") #"#000000",alpha=0)
cmapM.set_over("none")
cmapM.colorbar_extend="neither"
normm=BoundaryNorm(boundsM,cmapM.N) #len(bounds)-1)
############################################################
lon, lat, pox, poy = vs_locations_one("R_AMAZONAS_SOLIMOES_KM2637")
lon0=lon
lat0=lat
# flag= fflag[point]
# lon = flon[point]
# lat = flat[point]
# pox = fpox[point]
# poy = fpoy[point]
west, south = westsouth(lat,lon)
north = south + 10.0
east  = west + 10.0
cname0 = cname(lat,lon)
# get the dimesion of the map
dec=3
val=0.15
val1=0.10
val2=0.075
val3=0.05
val4=0.075
lllat = round_half_down(lat-val4,dec)
urlat = round_half_up(lat+val2,dec)
lllon = round_half_down(lon-val1,dec)
urlon = round_half_up(lon+val3,dec)
print (lllat, lllon, urlat, urlon)
# if abs(lllat-urlat) < val:
#     urlat=round_half_up(urlat+val,dec)
#     lllat=round_half_down(lllat-val,dec)
# if abs(lllon-urlon) < val:
#     urlon=round_half_up(urlon+val,dec)
#     lllon=round_half_down(lllon-val,dec)
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
# print (npix,":",spix,",",wpix,":",epix)
#=====================================
# high-resolution data
# print (cname0)
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
restag="3sec"
res=1.0/1200.0
visual=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".visual.bin"
# print (visual)
visual=np.fromfile(visual,np.int8).reshape(12000,12000)
#-----------------------------
# ax0 = fig.add_subplot(G[poy,pox])
# ax0.text(0.0,1.1,pname[point],va="center",ha="center",transform=ax0.transAxes,fontsize=14)
# flag_ch="Flag: %d"%(int(flag)) #lflag[point])
# ax0.text(0.5,0.05,flag_ch,va="center",ha="center",transform=ax0.transAxes,fontsize=6,color="w",zorder=113)
m = Basemap(projection='cyl',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon, lat_ts=0,resolution='c',ax=ax2)
try:
    # m.arcgisimage(service=maps[0], xpixels=1500, verbose=False)
    m.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=1200)
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
m.drawparallels([lllat,urlat], labels = [0,0,0,0], fontsize=4,linewidth=0,zorder=102)
m.drawmeridians([lllon,urlon], labels = [0,0,0,0], fontsize=4,linewidth=0,zorder=102)
data = ma.masked_less_equal(visual[npix:spix,wpix:epix],2)
im=m.imshow(data,interpolation="nearest",origin="upper",cmap=cmapM,norm=normm,zorder=110) # interpolation="nearest",origin="upper",
# print (lon,lat)
# m.scatter(lon,lat,s=0.5,marker="o",zorder=110,edgecolors="g", facecolors="g")#,transform=ccrs.PlateCarree()) #, 
# ax0.plot(lon ,lat ,color="g",marker="o",markersize=2,zorder=111) #fillstyle="none",
#================
orbit, lons, lats=sat_orbit(lllon,urlon,lllat,urlat)
for key in lons.keys(): #[0:10]:
    m.plot(lons[key],lats[key],color="brown",linewidth=5,alpha=0.5,zorder=110)
#================
# kx1= fkx1[point]
# ky1= fky1[point]
lat1 = south + 10.0 - res/2.0 - poy*res  
lon1 = west + res/2.0 + pox*res
# ax2.plot(lon1 ,lat1 ,color="y",marker="o",zorder=112)
ax2.plot(lon0 - 0.008 ,lat0 - 0.001,color="y",marker="o",zorder=112) #markersize=2,fillstyle="none",
ax2.plot(-68.491 ,-3.391 ,color="r",marker="o",zorder=113)
#================
plt.savefig("./fig/f02-satellite_schematic.png",dpi=800)
plt.savefig("./fig/f02-satellite_schematic.pdf",dpi=800)
plt.savefig("./fig/f02-satellite_schematic.jpg",dpi=800)