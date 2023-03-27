#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,Normalize,ListedColormap,BoundaryNorm
from matplotlib import colors
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import re
#=========================================
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
odir="/cluster/data6/menaka/AltiMaP/results"
if TAG=="HydroWeb":
    fname0=odir+"/HydroWeb/hydroweb_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/HydroWeb/hydroweb_cmf_linear_daily_wse_VIC_BC.nc"
    fname2=odir+"/HydroWeb/hydroweb_cmf_elediff_daily_wse_VIC_BC.nc"
    fname3=odir+"/HydroWeb/hydroweb_cmf_updown_daily_wse_VIC_BC.nc"
if TAG=="CGLS":
    fname0=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
    fname2=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
if TAG=="ICESat":
    fname0=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
    fname3=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
if TAG=="HydroSat":
    fname0=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
    fname2=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
############################################################
######################## original data #####################
nc0 = xr.open_dataset(fname0)
sfcelv_hydroweb0=nc0.sfcelv_hydroweb.values
sfcelv_cmf0=nc0.sfcelv_cmf.values
pname=nc0.name.values
lons=nc0.lon.values
lats=nc0.lat.values
basins=nc0.Basin.values
rivers=nc0.river.values
countries=nc0.country.values
sfcelv_hydroweb_max=nc0.sfcelv_hydroweb_max.values
sfcelv_hydroweb_min=nc0.sfcelv_hydroweb_min.values
sfcelv_cmf_max=nc0.sfcelv_cmf_max.values
sfcelv_cmf_min=nc0.sfcelv_cmf_min.values
disttomouth=nc0.disttomouth.values
elediff=nc0.elediff.values
flag=nc0.flag.values
nc0.close()
######################## interpolated data #####################
nc1 = xr.open_dataset(fname1)
pname1=nc1.name.values
sfcelv_cmf1=nc1.sfcelv_cmf.values
nc1.close()
######################## elevation differnce data #####################
nc2 = xr.open_dataset(fname2)
pname2=nc2.name.values
sfcelv_cmf2=nc2.sfcelv_cmf.values
nc2.close()
######################## upstream downstream data #####################
nc3 = xr.open_dataset(fname3)
pname3=nc3.name.values
sfcelv_cmf3=nc3.sfcelv_upstream_cmf.values
sfcelv_cmf4=nc3.sfcelv_downstream_cmf.values
nc3.close()
############################################################
# pnum=10 #len(pname)
#print np.shape(sfcelv_hydroweb)
colors=['xkcd:pastel blue','xkcd:teal','xkcd:aqua green','xkcd:dark pink','xkcd:purple','xkcd:magenta']
labels=["cmf oroginal","cmf interpolated","cmf ele diff",TAG]
#=============================
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
tag="3sec"
res=1.0/1200.0
#=============================
vmin=1.0
vmax=26.0
norm=Normalize(vmin=vmin,vmax=vmax)
bounds=np.arange(-0.5,26,1.0)
###
# sea=0
# land(undefined)=1
# land(defined in CaMa)=2 
# grid-box=3 
# catchment-boundary=5 
# channel=10 
# outlet-pixel=20 
# river-mouth=25
#
cmapL = matplotlib.colors.ListedColormap(['w','w','grey','k','w','k','w','y','w','w','blue','w','w','w','w','w','w','w','w','w','red', 'w','w','w','w','red'])
cmapL.set_under("none") #"#000000",alpha=0)
cmapL.set_over("none")
cmapL.colorbar_extend="neither"
norml=BoundaryNorm(bounds,cmapL.N) #len(bounds)-1)
############################################################
rivername0="CONGO" #"AMAZONAS"
stream0="CONGO" #"AMAZONAS"
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
kxlst=[]
kylst=[]
#-------------------------------------------
fname="./out/altimetry_"+mapname+"_test.txt"
# fname="./out/altimetry_"+mapname+"_20210518.txt"
#--
f=open(fname,"r")
lines=f.readlines()
for line in lines[1::]:
    line    = filter(None,re.split(" ",line))
    #print line
    num     = line[0]
    station = line[1]
    line2   = re.split("_",station)
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
    kx      = int(line[13])
    ky      = int(line[14])
    #-----------------------
    if riv != rivername0:
        continue
    if stream != stream0:
        continue
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
    kxlst.append(kx)
    kylst.append(ky)
#==========
pnum=len(pname)
#==========
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
#==========
if TAG=="HydroWeb":
    pdfname=odir+"/HydroWeb/hydroweb_cmf_hres_map.pdf"
if TAG=="CGLS":
    pdfname=odir+"/CGLS/cgls_cmf_hres_map.pdf"
if TAG=="ICESat":
    pdfname=odir+"/ICESat/icesat_cmf_hres_map.pdf"
if TAG=="HydroSat":
    pdfname=odir+"/HydroSat/hydrosat_cmf_hres_map.pdf"
#============================
with PdfPages(pdfname) as pdf:
    for point in np.arange(pnum):
        #prepare
        # cmf0=sfcelv_cmf0[:,point]
        # cmf1=sfcelv_cmf1[:,point]
        # cmf2=sfcelv_cmf2[:,point]
        # cmf3=sfcelv_cmf3[:,point]
        # cmf4=sfcelv_cmf4[:,point]
        # org=sfcelv_hydroweb0[:,point]
        # org0=ma.masked_where(sfcelv_hydroweb0[:,point]==-9999.0,sfcelv_hydroweb0[:,point])
        # locs=np.where(sfcelv_hydroweb0[:,point]!=-9999.0)[0]
        # org0=org0.compressed()
        ######################
        print ("=======================================")
        print (pname[point], lflag[point])
        hgt=11.69
        wdt=8.27
        fig=plt.figure(figsize=(wdt, hgt))
        #plt.title(pname[point][0],fontsize=12)
        G = gridspec.GridSpec(2,2)
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
        print (cname0)
        visual=CaMa_dir+"/map/"+mapname+"/"+tag+"/"+cname0+".visual.bin"
        # print (visual)
        visual=np.fromfile(visual,np.int8).reshape(12000,12000)
        #-----------------------------
        ax0 = fig.add_subplot(G[0,:])
        ax0.text(0.5,1.3,pname[point],va="center",ha="center",transform=ax0.transAxes,fontsize=14)
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
        m.drawparallels([lllat,urlat], labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
        m.drawmeridians([lllon,urlon], labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
        data = ma.masked_less_equal(visual[npix:spix,wpix:epix],2)
        im=m.imshow(data,interpolation="nearest",origin="upper",cmap=cmapL,norm=norml,zorder=110) # interpolation="nearest",origin="upper",
        print (lon,lat)
        # m.scatter(lon,lat,s=0.5,marker="o",zorder=110,edgecolors="g", facecolors="g")#,transform=ccrs.PlateCarree()) #, 
        ax0.plot(lon ,lat ,color="g",marker="o",markersize=7,zorder=111) #fillstyle="none",
        kx = kxlst[point]
        ky = kylst[point]
        lat0 = south + res/2.0 + 10.0 - ky*res  
        lon0 = west + res/2.0 + kx*res
        ax0.plot(lon0 ,lat0 ,color="r",marker="o",markersize=7,zorder=112) #fillstyle="none",
        print (kx,ky,lon0,lat0)
        #========================================================
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
    # set the file's metadata via the PdfPages object:
    d = pdf.infodict()
    d['Title'] = TAG+' high resolution map'
    d['Author'] = 'Menaka Revel'
    d['Subject'] = TAG+' high resolution map'
    d['Keywords'] = TAG+', high resolution, 3sec'
    d['CreationDate'] = datetime.datetime(2021, 5, 13)
    d['ModDate'] = datetime.datetime.today()