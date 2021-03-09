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

from read_patchMS import upstream
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
def upgrids(ix,iy,nums,nextxy,uparea,elevtn,nxtdst,rivseq,thr=50.0):
    # thr [m] threshold
    upXX=[ix]
    upYY=[iy]
    iXX=ix+1
    iYY=iy+1
    ele0=elevtn[iYY-1,iXX-1]
    ele1=elevtn[iYY-1,iXX-1]
    for _ in range(nums):
        nextX=nextxy[0]
        nextY=nextxy[1]
        ele0=elevtn[iYY-1,iXX-1]
        uXX, uYY = upstream(iXX,iYY,nextX.T,nextY.T,uparea.T)
        ele1=elevtn[uYY-1,uXX-1] 
        iXX=uXX
        iYY=uYY
        #---------
        if (ele1-ele0) > thr:
            break
        if iXX<0 or iYY<0:
            break
        if rivseq[iYY-1,iXX-1]<=1:
            break
        upXX.append(iXX-1)
        upYY.append(iYY-1)
    return list(reversed(upXX)), list(reversed(upYY))
#=============================
def dngrids(ix,iy,nums,nextxy,elevtn,thr=100.0):
    # thr [m] threshold
    dnXX=[]
    dnYY=[]
    iix=ix
    iiy=iy
    ele0=elevtn[iiy,iix]
    ele1=elevtn[iiy,iix]
    for _ in range(nums):
        iXX=nextxy[0,iiy,iix] - 1
        iYY=nextxy[1,iiy,iix] - 1
        ele1=elevtn[iiy,iix]
        ele0=elevtn[iYY,iXX]
        if ele1-ele0 > thr:
            break
        if iXX<0 or iYY<0:
            break
        iix=iXX
        iiy=iYY
        dnXX.append(iXX)
        dnYY.append(iYY)
    return dnXX, dnYY
#=============================
def river_reach(oxx,oyy,exx,eyy,nextxy,elevtn,nxtdst):
    dist=[0.0]
    elen=[elevtn[oyy,oxx]]
    ixx=oxx
    iyy=oyy
    d=nxtdst[oyy,oxx]*1e-3
    while (ixx!=exx and iyy!=eyy):
        ix=nextxy[0,iyy,ixx] - 1
        iy=nextxy[1,iyy,ixx] - 1
        ixx=ix
        iyy=iy
        if ixx<0 or iyy<0:
            break
        if ixx==exx and iyy==eyy:
            break
        elen.append(elevtn[iyy,ixx])
        dist.append(d)
        d=d+nxtdst[iyy,ixx]*1e-3
    dist.append(d)
    if nextxy[0,iyy,ixx]>0:
        elen.append(elevtn[iyy,ixx])
    else:
        elen.append(0.0)
    return dist, elen
#=============================
def reach_length(ixx,iyy,exx,eyy,nextxy,nxtdst):
    dist=0.0
    iix=ixx
    iiy=iyy
    d=nxtdst[iiy,iix]*1e-3
    while (iix!=exx and iiy!=eyy):
        ix=nextxy[0,iiy,iix] - 1
        iy=nextxy[1,iiy,iix] - 1
        iix=ix
        iiy=iy
        if iix<0 or iiy<0:
            break
        if iix==exx and iiy==eyy:
            break
        dist=d
        d=d+nxtdst[iiy,iix]*1e-3
    return dist
#=============================
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname="glb_01min"
syear=1992
eyear=2020
start=datetime.date(syear,1,1)
end=datetime.date(eyear,12,31)
days=(end-start).days + 1
linecolors = plt.cm.jet(np.linspace(0,1,10))
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
e08list={}
e96list={}
ediflist={}
#=============================
mapname=sys.argv[1]
CaMa_dir=sys.argv[2]
# mapname="glb_01min"
fname="./out/altimetry_"+mapname+".txt"
with open(fname, "r") as f:
    lines=f.readlines()
    for line in lines[1::]:
        line = filter(None,re.split(" ", line))
        #print line
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
            e08list[ixiy].append(egm08)
            e96list[ixiy].append(egm96)
            ediflist[ixiy].append(eled)
        else:
            pnames[ixiy]=[name]
            dataname[ixiy]=[dname]
            xxlist[ixiy]=[ix]
            yylist[ixiy]=[iy]
            lonlist[ixiy]=[lon]
            latlist[ixiy]=[lat]
            satlist[ixiy]=[sat.upper()]
            e08list[ixiy]=[egm08]
            e96list[ixiy]=[egm96]
            ediflist[ixiy]=[eled]
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
#=====================================
duplict=[]
closedict={}
upgrid={}
dngrid={}
nums=30
i=1
for ixiy in pnames.keys():
    if ixiy in duplict:
        continue
    duplict.append(ixiy)
    #print ixiy
    iflag="%05d"%(i)
    #print iflag
    closedict[iflag]=[]
    ix=xxlist[ixiy][0]
    iy=yylist[ixiy][0]
    uXX,uYY=upgrids(ix,iy,nums,nextxy,uparea,elevtn,nxtdst,rivseq)
    dXX,dYY=dngrids(ix,iy,nums,nextxy,elevtn)
    uXX.extend(dXX)
    uYY.extend(dYY)
    upgrid[iflag]=[uXX[0],uYY[0]]
    dngrid[iflag]=[uXX[-1],uYY[-1]]
    flag=0
    for j in range(len(uXX)):
        iXX=uXX[j]
        iYY=uYY[j]
        xxyy="%05d%05d"%(iXX,iYY)
        # if xxyy==ixiy:
        #     continue
        if xxyy in pnames.keys():
            #print xxyy
            closedict[iflag].append(xxyy)
            duplict.append(xxyy)
            flag=1
    if flag==1:
        i=i+1
#=====================================
# markers={"HydroWeb":"o","CGLS":"s","ICESat":"^","HydroSat":"X"}
# colors={"HydroWeb":"xkcd:reddy brown","CGLS":"xkcd:blue green","ICESat":"xkcd:pinkish","HydroSat":"xkcd:light urple"}
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
pdfname="./fig/WSE_observation_alongriver.pdf"
#=============================
with PdfPages(pdfname) as pdf:
    for iflag in closedict.keys():
        num=len(closedict[iflag])
        #print num
        if num<2:
            continue
        lons=[]
        lats=[]
        dnam=[]
        nams=[]
        lxx=[]
        lyy=[]
        le08=[]
        le96=[]
        ledf=[]
        for ixiy in closedict[iflag]:
            lons.extend(lonlist[ixiy])
            lats.extend(latlist[ixiy])
            dnam.extend(dataname[ixiy])
            nams.extend(pnames[ixiy])
            lxx.extend([int(ixiy[0:5])]*len(pnames[ixiy]))
            lyy.extend([int(ixiy[5::])]*len(pnames[ixiy]))
            le08.extend(e08list[ixiy])
            le96.extend(e96list[ixiy])
            ledf.extend(ediflist[ixiy])
        pnum=len(dnam)
        ###############################
        #         make figure         #
        ###############################
        hgt=11.69
        wdt=8.27
        linecolors = plt.cm.jet(np.linspace(0,1,pnum))
        fig=plt.figure(figsize=(wdt, hgt))
        #plt.title(pname[point][0],fontsize=12)
        G   = gridspec.GridSpec(3,2)
        ax0 = fig.add_subplot(G[0,0])  # google earth map
        ax1 = fig.add_subplot(G[1,0])  # longitudanal profile
        ax2 = fig.add_subplot(G[0:2,1])# Name list
        ax3 = fig.add_subplot(G[2,:])  # timeseries
        # get the dimesion of the map
        dec=2
        val=0.25
        lllat = round_half_down(min(lats)-val,dec)
        urlat = round_half_up(max(lats)+val,dec)
        lllon = round_half_down(min(lons)-val,dec)
        urlon = round_half_up(max(lons)+val,dec)
        #===============
        if abs(lllat-urlat) < val:
            urlat=round_half_up(urlat+val,dec)
            lllat=round_half_down(lllat-val,dec)
        if abs(lllon-urlon) < val:
            urlon=round_half_up(urlon+val,dec)
            lllon=round_half_down(lllon-val,dec)
        print (lllat, lllon, urlat, urlon)

        # if lllat==urlat:
        #     urlat=round_half_up(latlist[ixiy][0]+val,dec)
        #     lllat=round_half_down(latlist[ixiy][0]-val,dec)
        # if lllon==urlon:
        #     urlon=round_half_up(lonlist[ixiy][0]+val,dec)
        #     lllon=round_half_down(lonlist[ixiy][0]-val,dec)
        # dlat=urlat - lllat
        # dlon=urlon - urlon
        # if dlon>dlat:
        #     lllat=((lllat+urlat)/2.0)-dlon/2.0
        #     urlat=((lllat+urlat)/2.0)+dlon/2.0
        # else:
        #     lllon=((lllon+urlon)/2.0)-dlat/2.0
        #     urlon=((lllon+urlon)/2.0)+dlat/2.0
        # print lllat, lllon, urlat, urlon
        #------
        M = Basemap(resolution='h', projection='cyl',llcrnrlon=lllon, llcrnrlat=lllat, \
            urcrnrlon=urlon, urcrnrlat=urlat, ax=ax0)
        try:
            M.arcgisimage(service=maps[0], xpixels=1500, verbose=False)
        except:
            # Draw some map elements on the map
            M.drawcoastlines()
            M.drawstates()
            M.drawcountries()
            M.drawrivers(color='blue')
            M.fillcontinents(color="grey",lake_color='blue',zorder=99)
        M.drawparallels([lllat,urlat], labels = [1,0,0,0], fontsize=10,linewidth=0,zorder=102)
        M.drawmeridians([lllon,urlon], labels = [0,0,0,1], fontsize=10,linewidth=0,zorder=102)
        ##########
        # unilist=np.unique(np.array(dataname[ixiy]))
        # pnum=np.shape(np.unique(np.array(dataname[ixiy])))[0]
        # if pnum<2:
        #     continue
        obss=[]
        tags=[]
        lines=[]
        labels=[]
        j=0
        labeltag=[]
        oxx=upgrid[iflag][0]
        oyy=upgrid[iflag][1]
        exx=dngrid[iflag][0]
        eyy=dngrid[iflag][1]
        # oxx=lxx[0]
        # oyy=lyy[0]
        # exx=lxx[-1]
        # eyy=lyy[-1]
        dists,elevs=river_reach(oxx,oyy,exx,eyy,nextxy,elevtn,nxtdst)
        # print oxx, oyy, exx, eyy
        # print dists, elevs
        ax1.plot(dists,elevs,color="k")
        ax1.plot(dists,np.array(elevs)+15.0,color="k",linestyle="--")
        ax1.plot(dists,np.array(elevs)-10.0,color="k",linestyle="--")
        for i in np.arange(pnum):
            TAG=dnam[i] #dataname[ixiy][i]
            lon=lons[i]
            lat=lats[i]
            ixx=lxx[i]
            iyy=lyy[i]
            M.scatter(lon,lat,c=linecolors[i],s=50,marker=markers[TAG])
            #ax0.annotate(pname,xy=(lon,lat),fontsize=8,textcoords="offset points",xytext=(lon,lat),arrowprops=dict(arrowstyle="-"),zorder=111) #ha=ha,va=va,xycoords=transform,
            ax2.text(0.0,1.0-0.1*i,nams[i],va="center",ha="left",transform=ax2.transAxes,fontsize=10)
            locs,org = get_data(nams[i],TAG)
            dist=reach_length(oxx,oyy,ixx,iyy,nextxy,nxtdst)
            if TAG == "ICESat":
                print TAG, nams[i], np.mean(np.array(org)), "reach :", dist
                ax1.plot(dist,np.mean(np.array(org)),color=linecolors[i],marker=markers[TAG],markersize=5) #,fillstyle="none"
                ax1.plot(dist,np.mean(np.array(org))+ledf[i],color=linecolors[i],marker=markers[TAG],fillstyle="none",markersize=5)
            else:
                print TAG, nams[i], np.mean(np.array(org))+le08[i]-le96[i], "reach :", dist
                ax1.plot(dist,np.mean(np.array(org))+le08[i]-le96[i],color=linecolors[i],marker=markers[TAG],markersize=5) #,fillstyle="none"
                ax1.plot(dist,np.mean(np.array(org))+le08[i]-le96[i]+ledf[i],color=linecolors[i],marker=markers[TAG],fillstyle="none",markersize=5)
            lines.append(ax3.plot(locs,org,color=linecolors[i],label=nams[i],linestyle='-',linewidth=0.5)[0])
            labels.append(nams[i])
            #print locs,org
        #######################################
        #  print labels
        ax2.set_axis_off()
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax3.set_xlim(xmin=0,xmax=days+1)
        xxlab=np.arange(syear,eyear+1,2)
        dt=int(math.ceil(((eyear-syear)+2)/2.0))
        xxlist=np.linspace(0,days,dt,endpoint=True)
        ax3.set_xticks(xxlist)
        ax3.set_xticklabels(xxlab,fontsize=8)
        ax3.legend(lines,labels,ncol=1,loc='upper right',bbox_to_anchor=(1.1, 2.2))  
        #######################################
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