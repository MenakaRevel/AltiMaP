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
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes

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
mapname=sys.argv[1]
CaMa_dir=sys.argv[2]
ncpus=int(sys.argv[3])
# mapname="glb_01min"
fname="./out/altimetry_"+mapname+".txt"
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
#================================
#=== function for writing pdf ===
#================================
def write_pdf(inputlist):
    num   = inputlist[0]
    start = int(inputlist[1])
    last  = int(inputlist[2])
    pdfname="./fig/WSE_observation_river_network_"+num+".pdf"
    with PdfPages(pdfname) as pdf:
        for ixiy in pnames.keys():
            unilist=np.unique(np.array(dataname[ixiy]))
            pnum=np.shape(np.unique(np.array(dataname[ixiy])))[0]
            if pnum<2:
                continue
            hgt=11.69
            wdt=8.27
            fig=plt.figure(figsize=(wdt, hgt))
            #plt.title(pname[point][0],fontsize=12)
            G   = gridspec.GridSpec(3,2)
            ax00 = fig.add_subplot(G[0,0])
            ax01 = fig.add_subplot(G[0,1])
            ax1 = fig.add_subplot(G[1,:])
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
            print (lllat, lllon, urlat, urlon)
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
            obss=[]
            tags=[]
            lines=[]
            labels=[]
            j=0
            #print pnum
            z=0
            for i in np.arange(pnum):
                TAG=unilist[i] #dataname[ixiy][i]
                #print np.where(np.array(dataname[ixiy])==TAG)[0] #.index(TAG)
                repeatlist=np.where(np.array(dataname[ixiy])==TAG)[0]
                #print i, TAG, repeatlist
                obs_frg=[]
                k=0
                rflag=-9
                for j in repeatlist:
                    pname=pnames[ixiy][j]
                    lon=lonlist[ixiy][j]
                    lat=latlist[ixiy][j]
                    M.scatter(lon,lat,c=colors[TAG],s=20,marker=markers[TAG])
                    #ax0.annotate(pname,xy=(lon,lat),fontsize=8,textcoords="offset points",xytext=(lon,lat),arrowprops=dict(arrowstyle="-"),zorder=111) #ha=ha,va=va,xycoords=transform,
                    ax01.text(-0.1,1.0-0.1*z,pnames[ixiy][j],va="center",ha="left",transform=ax01.transAxes,fontsize=10)
                    z=z+1
                    #ax01.text(0.0,1.0-0.1*z,satlist[ixiy][j],va="center",ha="center",transform=ax0.transAxes,fontsize=14)
                    #z=z+1
                    #Observation
                    #print TAG, markers[TAG], colors[TAG]
                    locs,org = get_data(pname,TAG)
                    # try:
                    #     rflag=1
                    #     locs,org = get_data(pname,TAG)
                    # except:
                    #     rflag=0
                    #     continue
                    if k==0:
                        lines.append(ax1.plot(locs,org,color=colors[TAG],label=TAG,linestyle='-',linewidth=0.5,marker=markers[TAG],fillstyle="none",markersize=5)[0])
                    else:
                        ax1.plot(locs,org,color=colors[TAG],label=TAG,linestyle='-',linewidth=0.5,marker=markers[TAG],fillstyle="none",markersize=5)[0]
                    k=k+1
                obss.append(org)
                tags.append(TAG)
                labels.append(TAG)
                print (labels)
                # if rflag==1:
                #     obss.append(org)
                #     tags.append(TAG)
                #     labels.append(TAG)
                #     print labels
                # else:
                #     continue
            # ax00
            ax00.set_axis_off()
            ax00.spines['top'].set_visible(False)
            ax00.spines['right'].set_visible(False)
            ax00.spines['left'].set_visible(False)
            ax00.spines['bottom'].set_visible(False)
            #ax01
            ax01.set_axis_off()
            ax01.spines['top'].set_visible(False)
            ax01.spines['right'].set_visible(False)
            ax01.spines['left'].set_visible(False)
            ax01.spines['bottom'].set_visible(False)
            #ax1
            ax1.set_xlim(xmin=0,xmax=days+1)
            xxlab=np.arange(syear,eyear+1,5)
            dt=int(math.ceil(((eyear-syear)+2)/5.0))
            xxlist=np.linspace(0,days,dt,endpoint=True)
            ax1.set_xticks(xxlist)
            ax1.set_xticklabels(xxlab,fontsize=8)
            plt.legend(lines,labels,ncol=pnum+1,loc='upper right',bbox_to_anchor=(1.0, 1.13))  
            ######################################
            #boxplot
            ax2 = fig.add_subplot(G[2,0])
            flierprops = dict(marker='o', markerfacecolor='none', markersize=12,linestyle='none', markeredgecolor='k')
            boxprops = dict(color='grey')#facecolor='none'
            whiskerprops = dict(color='grey',linestyle="--")
            capprops = dict(color='grey')
            medianprops = dict(color='r')
            box=ax2.boxplot(obss,labels=labels,boxprops=boxprops,showfliers=False, \
                            whiskerprops=whiskerprops,capprops=capprops,medianprops=medianprops, \
                            notch=False, sym=None, vert=True, whis=1.5,positions=None, widths=None, \
                            patch_artist=True,bootstrap=None, usermedians=None, conf_intervals=None)#flierprops=flierprops,
            for patch, tag in zip(box['boxes'], labels):
                #patch.set_facecolor(colors[tag],alpha=0.5)
                patch.set(facecolor=colors[tag],alpha=0.5)
            ax2.set_ylabel('WSE $(m)$', color='k',fontsize=10)
            ax2.tick_params('y',labelsize=8, colors='k')
            ax2.tick_params('x',labelsize=8, colors='k')#,labelrotation=45)
            ax2.set_xticklabels(labels,rotation=0)
            ######################################
            #pdf
            ax3 = fig.add_subplot(G[2,1])
            #sns.distplot(obss[0], ax=ax3, hist=True, color="xkcd:cornflower", label="CaMa-Flood") #ax=ax3,
            for k,tag in enumerate(tags,start=0):
                sns.distplot(obss[k], ax=ax3, hist=True, color=colors[tag], label=tag) #ax=ax3,
            ax3.set_ylabel('density', color='k',fontsize=10)
            ax3.tick_params('y',labelsize=6, colors='k')
            ax3.set_xlabel('WSE $(m)$', color='k',fontsize=10)
            ax3.tick_params('x',labelsize=6, colors='k')
            #ax3.set_title("Histogram of Bias",fontsize=8)
            #ax3.set_xlim(xmin=-20.0,xmax=20.0)
            #ax3.text(0.01,0.95,"b",transform=ax3.transAxes,fontsize=8)
            ax3.legend(ncol=5, bbox_to_anchor=(0.0, -0.3), loc='lower center') #
            ######################################
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
            print ("============================")
        # set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Title'] = 'Comparison of HydroWeb, CGLS, HydroSat, ICESat altimetry data'
        d['Author'] = 'Menaka Revel'
        d['Subject'] = 'Comparison of altimetry observations'
        d['Keywords'] = 'HydroWeb, CGLS, HydroSat, ICESat, GRRATS'
        d['CreationDate'] = datetime.datetime(2021, 1, 25)
        d['ModDate'] = datetime.datetime.today()

#==============================================
pages=100.0
nums = int(len(pnames.keys())/pages) + 1
lastitem = len(pnames.keys())
inputlist=[]
for num in np.arange(0,nums):
    start = num*pages
    last  = (num + 1)*pages
    if last >= len(pnames.keys()):
        last = lastitem
    numch   = "%03d"%(num)
    startch = "%d"%(start)
    lastch  = "%d"%(last)
    inputlist.append([numch, startch, lastch])

#==========================
#== parallel writing pdf ==
#==========================
para_flag=1
# para_flag=0
#--
if para_flag==1:
    p=Pool(ncpus)
    p.map(write_pdf,inputlist)
    p.terminate()
else:
    map(write_pdf,inputlist)