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
import os
import errno
#=========================================
def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
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
def meanHydroWeb(station,egm08=0.0,egm96=0.0): 
    fname="/cluster/data6/menaka/HydroWeb/data/hydroprd_"+station+".txt"
    f=open(fname,"r")
    lines=f.readlines()
    f.close()
    head=33
    #--
    data=[] # WSE in [m]
    for line in lines[head::]:
        if line[0][0] == "#":
            continue
        line = re.split(" ",line)
        date = line[0]
        date = re.split("-",date)
        yyyy = int(date[0])
        mm   = int(date[1])
        dd   = int(date[2])
        wse  = float(line[2]) #+egm08-egm96
        data.append(wse)
    data=np.array(data)
    return np.mean(data), np.std(data)
#=============================
# sfcelv
syear=2002
eyear=2013
#=========================================
#TAG="CGLS"
TAG="HydroWeb"
# TAG="ICESat"
# TAG="HydroSat"
# # # #=========================================
# # # # odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
# # # # fname = odir+"/hydroweb_cmf_daily_wse_VIC_BC.nc"
# # # odir="/cluster/data6/menaka/Altimetry/results"
# # # if TAG=="HydroWeb":
# # #     fname0=odir+"/HydroWeb/hydroweb_cmf_daily_wse_VIC_BC.nc"
# # #     fname1=odir+"/HydroWeb/hydroweb_cmf_linear_daily_wse_VIC_BC.nc"
# # #     fname2=odir+"/HydroWeb/hydroweb_cmf_elediff_daily_wse_VIC_BC.nc"
# # #     fname3=odir+"/HydroWeb/hydroweb_cmf_updown_daily_wse_VIC_BC.nc"
# # # if TAG=="CGLS":
# # #     fname0=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
# # #     fname1=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
# # #     fname2=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
# # # if TAG=="ICESat":
# # #     fname0=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
# # #     fname1=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
# # #     fname3=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
# # # if TAG=="HydroSat":
# # #     fname0=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
# # #     fname1=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
# # #     fname2=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
# # # ############################################################
# # # ######################## original data #####################
# # # nc0 = xr.open_dataset(fname0)
# # # sfcelv_hydroweb0=nc0.sfcelv_hydroweb.values
# # # sfcelv_cmf0=nc0.sfcelv_cmf.values
# # # pname=nc0.name.values
# # # lons=nc0.lon.values
# # # lats=nc0.lat.values
# # # basins=nc0.Basin.values
# # # rivers=nc0.river.values
# # # countries=nc0.country.values
# # # sfcelv_hydroweb_max=nc0.sfcelv_hydroweb_max.values
# # # sfcelv_hydroweb_min=nc0.sfcelv_hydroweb_min.values
# # # sfcelv_cmf_max=nc0.sfcelv_cmf_max.values
# # # sfcelv_cmf_min=nc0.sfcelv_cmf_min.values
# # # disttomouth=nc0.disttomouth.values
# # # elediff=nc0.elediff.values
# # # flag=nc0.flag.values
# # # nc0.close()
# # # ######################## interpolated data #####################
# # # nc1 = xr.open_dataset(fname1)
# # # pname1=nc1.name.values
# # # sfcelv_cmf1=nc1.sfcelv_cmf.values
# # # nc1.close()
# # # ######################## elevation differnce data #####################
# # # nc2 = xr.open_dataset(fname2)
# # # pname2=nc2.name.values
# # # sfcelv_cmf2=nc2.sfcelv_cmf.values
# # # nc2.close()
# # # ######################## upstream downstream data #####################
# # # nc3 = xr.open_dataset(fname3)
# # # pname3=nc3.name.values
# # # sfcelv_cmf3=nc3.sfcelv_upstream_cmf.values
# # # sfcelv_cmf4=nc3.sfcelv_downstream_cmf.values
# # # nc3.close()
# # # ############################################################
mkdir("./fig/"+TAG)
mkdir("./fig/"+TAG+"/along_river_VS")
# pnum=10 #len(pname)
#print np.shape(sfcelv_hydroweb)
colors=['xkcd:pastel blue','xkcd:teal','xkcd:aqua green','xkcd:dark pink','xkcd:purple','xkcd:magenta']
labels=["cmf oroginal","cmf interpolated","cmf ele diff",TAG]
#=============================
mapname="glb_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
tag="3sec"
res=1.0/1200.0
upthr=10.0
dwthr=10.0
#=============================
# Read the CMF variables
if mapname == 'glb_15min':
    nx      = 1440
    ny      = 720
    ny_     = 640
elif mapname == 'glb_06min':
    nx      = 3600
    ny      = 1800
    ny_     = 1500
elif mapname == 'glb_01min':
    nx      = 21600
    ny      = 10800
    ny_     = 10800
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
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
# rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
rivseq = np.fromfile(rivseq,np.int32).reshape(ny,nx)
#---
nextX=nextxy[0]
nextY=nextxy[1]
#=============================
meanWSE_VICBC="/cluster/data6/menaka/ensemble_org/CaMa_out/GLBVIC_BC001/sfcelv_mean2000-2013.bin"
meanWSE_VICBC=np.fromfile(meanWSE_VICBC,np.float32).reshape(ny_,nx)
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
rivernames=["CONGO"] #["AMAZONAS","CONGO","AMAZONAS","LENA","GANGES-BRAHMAPUTRA","MEKONG","MISSISSIPPI","IRRAWADDY"]
streams=["CONGO"] #["XINGU","CONGO","AMAZONAS","LENA","BRAHMAPUTRA","MEKONG","MISSISSIPPI","IRRAWADDY"]
for j in np.arange(len(rivernames)):
    rivername0 = rivernames[j]
    stream0 = streams[j]
    ############################################################
    nums=[]
    river=[]
    pname=[]
    rvlen=[]
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
    # fname="./out/altimetry_"+mapname+"_test.txt"
    # fname="./out/altimetry_"+mapname+"_20210602.txt"
    # fname="./out/altimetry_"+mapname+"_20210531.txt"
    fname="./out/altimetry_"+mapname+"_20210617.txt"
    # fname="./tmp.txt"
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
        # print (line2[-1])
        rlen    = float(line2[-1][2::])
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
        #-----------------------
        # print (riv, stream, rlen)
        nums.append(num)
        river.append(riv)
        pname.append(station)
        rvlen.append(rlen)
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
    #-------------------------------------------
    rvlen=np.array(rvlen)
    pname=np.array(pname)
    ylist=np.array(ylist)
    xlist=np.array(xlist)
    egm08=np.array(egm08)
    egm96=np.array(egm96)
    #--
    order=np.argsort(rvlen)[::-1]
    #--
    rvlen00=np.amax(rvlen) - rvlen[order]
    pname00=pname[order]
    ylist00=ylist[order]
    xlist00=xlist[order]
    egm0800=egm08[order]
    egm9600=egm96[order]
    pnum=len(pname00)
    #----------------
    mean_elevation=[]
    mean_WSE=[]
    std_WSE=[]
    cmf_WSE=[]
    riv_hgt=[]
    dist=0.0
    for point in np.arange(pnum):
        # print (pname00[point], rvlen00[point], meanHydroWeb(pname00[point],egm96=egm9600[point],egm08=egm0800[point]), elevtn[ylist00[point],xlist00[point]], 
        #         rivseq[ylist00[point],xlist00[point]], meanWSE_VICBC[ylist00[point],xlist00[point]])
        meanW, stdW = meanHydroWeb(pname00[point],egm96=egm9600[point],egm08=egm0800[point])
        mean_WSE.append(meanW)
        std_WSE.append(stdW)
        mean_elevation.append(elevtn[ylist00[point],xlist00[point]])
        cmf_WSE.append(meanWSE_VICBC[ylist00[point],xlist00[point]])
        riv_hgt.append(rivhgt[ylist00[point],xlist00[point]])
    mean_WSE=np.array(mean_WSE)
    std_WSE=np.array(std_WSE)
    cmf_WSE=np.array(cmf_WSE)
    mean_elevation=np.array(mean_elevation)
    riv_hgt=np.array(riv_hgt)
    #-------------------------------------------
    hgt=11.69*(1.0/2.0)
    wdt=8.27
    fig=plt.figure(figsize=(wdt, hgt))
    #plt.title(pname[point][0],fontsize=12)
    G  = gridspec.GridSpec(1,1)
    ax = fig.add_subplot(G[0,0])
    ax.set_title(rivername0+"-"+stream0,fontsize=14)
    # ax.plot(rvlen00,mean_WSE,color="k",label=TAG,linestyle='none',linewidth=0,marker="o",fillstyle="none",markersize=5)
    ax.errorbar(rvlen00,mean_WSE,yerr=std_WSE,color="k",label=TAG,linestyle='none',linewidth=0.5,marker="o",fillstyle="none",markersize=2,
                ecolor="k",elinewidth=0.5,capsize=1.0)
    ax.plot(rvlen00,mean_elevation,color="g",label="CaMa-Flood")
    ax.plot(rvlen00,mean_elevation + 10.0,color="xkcd:bluish green",label="CaMa-Flood",linestyle='--',linewidth=0.5)
    ax.plot(rvlen00,mean_elevation - 10.0,color="xkcd:bluish green",label="CaMa-Flood",linestyle='--',linewidth=0.5)
    # ax.plot(rvlen00,mean_elevation - riv_hgt,color="purple",label="CaMa-Flood",linestyle='--',linewidth=0.5) #"xkcd:bluish green"
    ax.plot(rvlen00,cmf_WSE,color="r",label="CaMa-Flood",linestyle='-',linewidth=0.5) #xkcd:bright blue
    # plot the unacceptble VS
    slp_thr=1.0*abs((mean_elevation[0]-mean_elevation[-1])/rvlen00[-1])
    print slp_thr
    slope0=0.0
    slope_org = [0.0] #[abs((mean_elevation[i]-mean_elevation[i-1])/(rvlen00[i-1]-rvlen00[i])) for i in range(1,pnum)]
    print slope_org 
    slopes=[]
    slopes0=[]
    slp_len=[]
    for point in range(pnum):
        if point < 1:
            continue
        #-------------------
        uppoint1 = point - 1
        dwpoint1 = point + 1
        if mean_elevation[point] + upthr < mean_WSE[point] or mean_elevation[point] - dwthr > mean_WSE[point]:
            ax.plot(rvlen00[point],mean_WSE[point],color="r",linestyle='none',linewidth=0.5,marker="o",fillstyle="none",markersize=10)
        #-------
        slope_obs = abs((mean_WSE[uppoint1]-mean_WSE[point])/(rvlen00[point]-rvlen00[uppoint1]))
        slope_cmf = abs((mean_elevation[uppoint1]-mean_elevation[point])/(rvlen00[point]-rvlen00[uppoint1]))
        slopes.append(slope_cmf)
        slopes.append(slope_cmf)
        slope_arr=np.array(slopes)
        mean_slope=np.mean(slope_arr)
        slope_org.append(mean_slope)
        if abs(mean_slope - slope_cmf) > slp_thr:
            print "="*50
            slopes=[]
            slp_len.append(rvlen00[point])
        # print "%69s\t%12.5e\t%12.5e\t%12.5e"%(pname00[point], slope_obs, slope_cmf, slope_cmf-slope0)
        slope0=slope_cmf
        
        # if mean_elevation[point] + upthr < mean_WSE[point] or mean_elevation[point] - dwthr > mean_WSE[point]:
        #     print (pname00[point],rivername0,rvlen00[point],mean_WSE[point]) #,mean_WSE[point],mean_WSE[dwpoint1])
        #     ax.plot(rvlen00[point],mean_WSE[point],color="r",linestyle='none',linewidth=0.5,marker="o",fillstyle="none",markersize=10)
        
        # if point < 2:
        #     continue
        # if point > pnum -2:
        #     continue
        # #-------------------
        # uppoint1 = point - 1
        # dwpoint1 = point + 1
        # uppoint2 = point - 2
        # dwpoint2 = point + 2
        # flag = 0
        # # print (mean_WSE[uppoint],mean_WSE[point],mean_WSE[dwpoint])
        # if mean_WSE[uppoint1] + upthr < mean_WSE[point] or mean_WSE[dwpoint1] - dwthr > mean_WSE[point]:
        #     if mean_WSE[uppoint2] < mean_WSE[point] or mean_WSE[dwpoint2] > mean_WSE[point]:
        #         print (pname00[point],rivername0,rvlen00[point],mean_WSE[uppoint1],mean_WSE[point],mean_WSE[dwpoint1])
        #         ax.plot(rvlen00[point],mean_WSE[point],color="r",linestyle='none',linewidth=0.5,marker="o",fillstyle="none",markersize=10)
    #--
    slp00=[abs(slope_org[i]-slope_org[i-1]) for i in range(1,pnum)]
    slp_len=(slp00 > slp_thr)*rvlen00[1::]
    # slp_len=np.array(slp_len)
    # ax.vlines(x=slp_len,ymin=0.0,ymax=max(mean_elevation + 10.0), colors="yellow")
    ax.set_ylabel('Mean Elevation $(m)$', color='k',fontsize=10)
    ax.tick_params('y',labelsize=6, colors='k')
    ax.set_xlabel('Distance along river\nto river mouth $(km)$', color='k',fontsize=10)
    ax.tick_params('x',labelsize=6, colors='k')
    ax.set_xlim(xmin=0.0) #,xmax=1500.0)
    # ax.set_ylim(ymin=250.0)
    print (rivername0+"-"+stream0)
    plt.savefig("./fig/"+TAG+"/along_river_VS/"+rivername0+"-"+stream0+".png",dpi=500) 