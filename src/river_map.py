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
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import matplotlib.patches as mpatches
import re
import os
import sys
import errno

import read_hydroweb as hweb
import read_cgls as cgls
import read_hydrosat as hsat
import read_icesat as isat
import read_grrats as grt
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
def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5) / multiplier
#=============================
def round_half_down(n, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(n*multiplier - 0.5) / multiplier
#=============================
def vec_par(LEVEL,ax=None):
    sup=2
    w=0.02
    width=0.5
    ax=ax or plt.gca()
    txt="tmp_%02d.txt"%(LEVEL)
    os.system("./src/print_rivvec tmp_00.txt 1 "+str(LEVEL)+" > "+txt)
    width=(float(LEVEL)**sup)*w
    with open(txt,"r") as f:
        lines = f.readlines()
    #---
    for line in lines:
        line    = filter(None, re.split(" ",line))
        lon1 = float(line[0])
        lat1 = float(line[1])
        lon2 = float(line[3])
        lat2 = float(line[4])

        if lon1-lon2 > 180.0:
            print lon1, lon2
            lon2=180.0
        elif lon2-lon1> 180.0:
            print lon1,lon2
            lon2=-180.0

        colorVal="xkcd:azure"
        # print lon1,lon2,lat1,lat2
        plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=ax00)
#=============================
def plot_ax(lon1,lon2,lat1,lat2,width,colorVal,ax=None,alpha=1):
    ax=ax or plt.gca()
    return ax.plot([lon1,lon2],[lat1,lat2],color=colorVal,linewidth=width,zorder=105,alpha=alpha) #transform=ccrs.PlateCarree(),
#=============================
syear=1992
eyear=2020
start=datetime.date(syear,1,1)
end=datetime.date(eyear,12,31)
days=(end-start).days + 1
#=============================
idnames={}
pnames={}
dataname={}
listnum={}
xxlist={}
yylist={}
satlist={}
lonlist={}
latlist={}
#=============================
# mapname="glb_01min"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# print sys.argv
mapname=sys.argv[1]
CaMa_dir=sys.argv[2]
fname="./out/altimetry_"+mapname+"_new.txt"
with open(fname, "r") as f:
    lines=f.readlines()
    for line in lines[1::]:
        line = filter(None,re.split(" ", line))
        # print line
        Id   = line[0].strip()
        name = line[1].strip()
        dname= line[2].strip()
        lon  = float(line[3])
        lat  = float(line[4])
        ix   = int(line[5]) 
        iy   = int(line[6])
        eled = float(line[7])
        egm08= float(line[8])
        egm96= float(line[9])
        sat  = line[10].strip()
        ixiy="%05d%05d"%(ix,iy)
        #print ixiy
        if ixiy in pnames.keys():
            pnames[ixiy].append(name)
            dataname[ixiy].append(dname)
            xxlist[ixiy].append(ix)
            yylist[ixiy].append(iy)
            lonlist[ixiy].append(lon)
            latlist[ixiy].append(lat)
            satlist[ixiy].append(sat.upper())
        else:
            pnames[ixiy]=[name]
            dataname[ixiy]=[dname]
            xxlist[ixiy]=[ix]
            yylist[ixiy]=[iy]
            lonlist[ixiy]=[lon]
            latlist[ixiy]=[lat]
            satlist[ixiy]=[sat.upper()]
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
mkdir("./fig")
pdfname="./fig/WSE_observation_river_network.pdf"
#=============================
with PdfPages(pdfname) as pdf:
    for ixiy in pnames.keys():
        unilist=np.unique(np.array(dataname[ixiy]))
        pnum=np.shape(np.unique(np.array(dataname[ixiy])))[0]
        # if pnum<2:
        #     continue
        hgt=11.69
        wdt=8.27
        fig=plt.figure(figsize=(wdt, hgt))
        #plt.title(pname[point][0],fontsize=12)
        G   = gridspec.GridSpec(3,2)
        ax00 = fig.add_subplot(G[0:2,0::])
        #ax01 = fig.add_subplot(G[0,1])
        ax11 = fig.add_subplot(G[2,:])
        # get the dimesion of the map
        dec=2
        val=0.25
        lllat = round_half_down(min(latlist[ixiy])-val,dec)
        urlat = round_half_up(max(latlist[ixiy])+val,dec)
        lllon = round_half_down(min(lonlist[ixiy])-val,dec)
        urlon = round_half_up(max(lonlist[ixiy])+val,dec)
        if abs(lllat-urlat) < val:
            urlat=round_half_up(urlat+val,dec)
            lllat=round_half_down(lllat-val,dec)
        if abs(lllon-urlon) < val:
            urlon=round_half_up(urlon+val,dec)
            lllon=round_half_down(lllon-val,dec)
        print lllat, lllon, urlat, urlon
        M = Basemap(resolution='h', projection='cyl',llcrnrlon=lllon, llcrnrlat=lllat, \
            urcrnrlon=urlon, urcrnrlat=urlat, ax=ax00)
        try:
            M.arcgisimage(service=maps[0], xpixels=1500, verbose=False)
        except:
            # Draw some map elements on the map
            M.drawcoastlines()
            M.drawstates()
            M.drawcountries()
            M.drawrivers(color='blue')
        M.drawparallels([lllat,urlat], labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
        M.drawmeridians([lllon,urlon], labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
        #####
        box="%f %f %f %f"%(lllon,urlon,urlat,lllat) 
        # print box
        os.system("./src/txt_vector "+box+" "+CaMa_dir+" "+mapname+" > tmp_00.txt") 
        map(vec_par,np.arange(1,10+1,1))
        os.system("rm -r tmp*.txt")
        for i in np.arange(pnum):
            TAG=unilist[i] 
            repeatlist=np.where(np.array(dataname[ixiy])==TAG)[0]
            for j in repeatlist:
                pname=pnames[ixiy][j]
                lon=lonlist[ixiy][j]
                lat=latlist[ixiy][j]
                M.scatter(lon,lat,c=colors[TAG],s=20,marker=markers[TAG],zorder=110)
                ax11.text(-0.1,1.1-0.1*j,pnames[ixiy][j],va="center",ha="left",transform=ax11.transAxes,fontsize=10)
        ax11.set_axis_off()
        ax11.spines['top'].set_visible(False)
        ax11.spines['right'].set_visible(False)
        ax11.spines['left'].set_visible(False)
        ax11.spines['bottom'].set_visible(False)
        ######################################
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
        print "============================"
    # set the file's metadata via the PdfPages object:
    d = pdf.infodict()
    d['Title'] = 'Comparison of HydroWeb, CGLS, HydroSat, ICESat along the river'
    d['Author'] = 'Menaka Revel'
    d['Subject'] = 'Comparison of altimetry observations'
    d['Keywords'] = 'HydroWeb, CGLS, HydroSat, ICESat, GRRATS'
    d['CreationDate'] = datetime.datetime(2021, 1, 25)
    d['ModDate'] = datetime.datetime.today()