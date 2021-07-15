#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import datetime
import numpy as np
from numpy import ma
from numpy import dot
from numpy.linalg import solve
from numpy.polynomial.polynomial import Polynomial as P, polyvander as V
import re
import sys
import os
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import warnings;warnings.filterwarnings('ignore')
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes

sys.path.append("./src")
from read_patchMS import upstream
from river_function import river_profile
import read_hydroweb as hweb
import read_cgls as cgls
import read_hydrosat as hsat
import read_icesat as isat
import read_grrats as grt
#======================================================================
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
#==================
def clsq(A, b, C, d, M= 1e5):
    """A simple constrained least squared solution of Ax= b, s.t. Cx= d,
    based on the idea of weighting constraints with a largish number M."""
    return solve(dot(A.T, A)+ M* dot(C.T, C), dot(A.T, b)+ M* dot(C.T, d))
#==================
def cpf(x, y, x_c, y_c, n, M= 1e5):
    """Constrained polynomial fit based on clsq solution."""
    return P(clsq(V(x, n), y, V(x_c, n), y_c, M))
#==================
# calculate R^2
def rsq(y1,y2):
    yresid= y1 - y2
    SSresid = np.sum(yresid**2)
    SStotal = len(y1) * np.var(y1)
    r2 = 1 - SSresid/SStotal
    return r2
#=============================
def river_along_profile(ix,iy,west,south,csize,nx,ny,hiresmap):
    length,elevation,k = river_profile(ix,iy,west,south,csize,nx,ny,hiresmap)
    # print (length, elevation, k)
    return length[0:k] , elevation[0:k]
#=========================================
def westsouth(lat,lon):
    return float(int(math.floor(lon/10.0)*10)), float(int(math.floor(lat/10.0)*10))
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
#==================
def polyfit_with_fixed_points(n, x, y, xf, yf) :
    mat = np.empty((n + 1 + len(xf),) * 2)
    vec = np.empty((n + 1 + len(xf),))
    x_n = x**np.arange(2 * n + 1)[:, None]
    yx_n = np.sum(x_n[:n + 1] * y, axis=1)
    x_n = np.sum(x_n, axis=1)
    idx = np.arange(n + 1) + np.arange(n + 1)[:, None]
    mat[:n + 1, :n + 1] = np.take(x_n, idx)
    xf_n = xf**np.arange(n + 1)[:, None]
    mat[:n + 1, n + 1:] = xf_n / 2
    mat[n + 1:, :n + 1] = xf_n.T
    mat[n + 1:, n + 1:] = 0
    vec[:n + 1] = yx_n
    vec[n + 1:] = yf
    params = np.linalg.solve(mat, vec)
    return params[:n + 1]
#==================
def AIC_val(y,yfit,d):
    resid = sum(np.sqrt((y-yfit)**2))
    n=float(len(y))
    p=float(d+1)
    #--
    AIC=n*math.log(resid/n)+2*p
    # AIC=-2.0*math.log(resid/n)+2.0*p
    return AIC
#==================
colors=["xkcd:sea blue","xkcd:tangerine","xkcd:dark pink","xkcd:pinkish","xkcd:light urple"]
#======================================================================
rivername0="MEKONG" #"LENA" #"CONGO" #"AMAZONAS" #"AMAZONAS"
stream0="MEKONG" #"LENA" #"CONGO" #"AMAZONAS" #
# dataname="HydroWeb"
odir="/cluster/data6/menaka/Altimetry/fig/river_profile"
TAG="HydroWeb"
mapname="glb_06min"
restag="3sec"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
#======================================================================
os.system("mkdir -p "+odir)
#======================================================================
res=1.0/1200.0
nx =12000
ny =12000
hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
if restag == "3sec":
    res=1.0/1200.0
    nx =12000
    ny =12000
    hiresmap=CaMa_dir+"/map/"+mapname+"/"+restag+"/"
#=============================
# Read the CMF variables
if mapname == 'glb_15min':
    nXX      = 1440
    nYY      = 720
elif mapname == 'glb_06min':
    nXX      = 3600
    nYY      = 1800
elif mapname == 'glb_01min':
    nXX      = 21600
    nYY      = 10800
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
nextxy = np.fromfile(nextxy,np.int32).reshape(2,nYY,nXX)
# rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
# rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
# rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(nYY,nXX)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,nYY,nXX)
uparea = np.fromfile(uparea,np.float32).reshape(nYY,nXX)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(nYY,nXX)
rivseq = np.fromfile(rivseq,np.int32).reshape(nYY,nXX)
#======================================================================
obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"+mapname+"_20210709.txt"
fname=obstxt
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
############################################################
with open(fname,"r") as f:
    lines=f.readlines()
#==========================================
def mk_pdf(inputlist):
    # rivername0=rivernames[nriver]
    # stream0=streams[nriver]
    rivername0=inputlist[0]
    stream0=inputlist[1]
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
        # print (riv,station,kx,ky)
        if riv != rivername0:
            continue
        if riv == "AMAZONAS":
            if stream != stream0:
                continue
        elif riv == "CONGO":
            if stream != stream0:
                continue
        else:
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
        # print (riv,station)
    #=============================
    pnum=len(pname)
    #=====================================
    # #=============================
    # if TAG=="HydroWeb":
    #     pdfname=odir+"/HydroWeb/high_res/hydroweb_cmf_hres_map_"+rivername0+"-"+stream0+".pdf"
    # if TAG=="CGLS":
    #     pdfname=odir+"/CGLS/high_res/cgls_cmf_hres_map_"+rivername0+"-"+stream0+".pdf"
    # if TAG=="ICESat":
    #     pdfname=odir+"/ICESat/high_res/icesat_cmf_hres_map_"+rivername0+"-"+stream0+".pdf"
    # if TAG=="HydroSat":
    #     pdfname=odir+"/HydroSat/high_res/hydrosat_cmf_hres_map_"+rivername0+"-"+stream0+".pdf"
    pdfname=odir+"/river_prifles_"+rivername0+".pdf"
    #============================
    #############################
    with PdfPages(pdfname) as pdf:
        for point in np.arange(0,pnum):
            print ("=======================================")
            print (point, pname[point], lflag[point])
            hgt=11.69
            wdt=8.27
            fig=plt.figure(figsize=(wdt, hgt))
            #plt.title(pname[point][0],fontsize=12)
            G = gridspec.GridSpec(3,2)
            # lon = lons[point]
            # lat = lats[point]
            lon = lons[point]
            lat = lats[point]
            west, south = westsouth(lat,lon)
            north = south + 10.0
            east  = west + 10.0
            cname0 = cname(lat,lon)
            #========================================================
            # high-resolution data
            # print (cname0)
            visual=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".visual.bin"
            # print (visual)
            visual=np.fromfile(visual,np.int8).reshape(12000,12000)
            ele1m=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".elevtn.bin"
            # print (ele1m)
            ele1m=np.fromfile(ele1m,np.float32).reshape(12000,12000)
            #-----------------------------
            ax0 = fig.add_subplot(G[0,0])
            ax0.text(1.0,1.3,pname[point],va="center",ha="center",transform=ax0.transAxes,fontsize=14)
            flag_ch="flag: %d"%(lflag[point])
            ax0.text(1.0,1.1,flag_ch,va="center",ha="center",transform=ax0.transAxes,fontsize=14)
            ax0.axis("off")
            #----------------------------
            try:
                length=[1.0]
                elevation=[1.0]
                kx = kxlst[point]
                ky = kylst[point]
                ix = xlist[point]
                iy = ylist[point]
                print ("calculate river profile...........",kx,ky,visual[ky-1,kx-1])
                ax1 = fig.add_subplot(G[1,0])
                # print (kx,ky,west,south,res,nx,ny)
                length, elevation = river_along_profile(kx,ky,west,south,res,nx,ny,hiresmap)
                # print (length, elevation)
                ax1.plot(length,elevation,color="grey",linestyle='-',linewidth=3.0)
                disttom=ldtom[point]
                ax1.plot([length[-1]-disttom*1e3],[ele1m[ky-1,kx-1]],color="r",marker="o",markersize=7,linestyle='none',linewidth=0.0)
                ax1.plot([0,length[-1]],[elevation[0],elevation[-1]],color="g",linestyle='--',linewidth=0.5)
                elv_char="elevation diff=%5.3f"%(ele1m[ky-1,kx-1]-elevtn[iy,ix])
                ax1.text(1.4,0.0,elv_char,va="center",ha="center",color="k",transform=ax1.transAxes,fontsize=14)
                #-----------
                xf = np.array([length[0],length[-1]])
                yf = np.array([elevation[0],elevation[-1]])
                M  = 1e5
                lrsq=[]
                lAIC=[]
                print ("Fitting the polynomial..........")
                for n in np.arange(1,5+1):
                    try:
                        params = cpf(length, elevation, xf, yf, n, M)
                        # print (params)
                        ax1.plot(length,params(length),color=colors[n-1],linestyle='--',linewidth=0.5)
                        r_squared = rsq(elevation,params(length))
                        AICval = AIC_val(elevation,params(length),n)
                        lrsq.append(r_squared)
                        lAIC.append(AICval)
                        rsq_char="$r^2$ : %3.2f"%(r_squared)
                        AIC_char="$A\^IC$ : %3.2f"%(AICval)
                        print ("Fitting the model.........",n, r_squared, AICval)
                        ax1.text(1.2,1.1-n*0.1,rsq_char,va="center",ha="center",color=colors[n-1],transform=ax1.transAxes,fontsize=14)
                        ax1.text(1.7,1.1-n*0.1,AIC_char,va="center",ha="center",color=colors[n-1],transform=ax1.transAxes,fontsize=14)
                    except:
                        lrsq.append(-9999)
                        lAIC.append(1e20)
                        print ("cannot fit the model")
                #==================
                print ("Calculating statistics......")
                lrsq=np.array(lrsq)
                lAIC=np.array(lAIC)
                # replace nan to number
                print ("Removing nan values......")
                lrsq=np.nan_to_num(lrsq)#, nan=-9999, posinf=-9999, neginf=-9999)#, nan=-9999)
                lAIC=np.nan_to_num(lAIC)#, nan=1e20, posinf=1e20, neginf=1e20)
                lrsq=ma.masked_equal(lrsq,0.0).filled(-9999)
                lAIC=ma.masked_equal(lAIC,0.0).filled(1e20)
                #------------------
                print ("Calculation min values......", lrsq, lAIC)
                AIC_min=np.min(lAIC)
                AIC_mloc=np.argmin(lAIC)
                print (AIC_min, AIC_mloc)
                # compare with linear treand
                deltaAIC=lAIC[0]-AIC_min
                dAIC_char="$\{delta}A\^IC : %3.2f"%(deltaAIC)
                AICm_char="$A\^IC_min : %5.2f"%(AIC_min)
                mdle_char="Model Order: %d"%(AIC_mloc+1)
                print ("+++Statistics++++\n"+dAIC_char+"\n"+AICm_char+"\n"+mdle_char)
                ax1.text(1.5,-0.5,dAIC_char,va="center",ha="center",color="k",transform=ax1.transAxes,fontsize=14)
                ax1.text(1.5,-0.3,AICm_char,va="center",ha="center",color="k",transform=ax1.transAxes,fontsize=14)
                ax1.text(1.5,-0.1,mdle_char,va="center",ha="center",color="k",transform=ax1.transAxes,fontsize=14)
            except:
                length=[0.0]
                elevation=[0.0]
            # # print (length, elevation)
            if visual[ky-1,kx-1]!=10: # or visual[ky-1,kx-1]!=20:
                length=[0.0]
                elevation=[0.0]
            #========================================================
            locs,org = get_data(pname[point],TAG,egm08=egm08[point],egm96=egm96[point])
            ax1.axhline(y=np.mean(org),xmin=0.0,xmax=length[-1],color="blue",linestyle='--',linewidth=1.0)
            ax1.set_xlim(xmin=0,xmax=length[-1])
            ax1.set_xlabel("Distance $(m)$")
            ax1.set_ylabel("Elevation $(m EGM96)$")
            #========================================================
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
        # set the file's metadata via the PdfPages object:
        d = pdf.infodict()
        d['Title'] = ' high resolution river profile '+rivername0
        d['Author'] = 'Menaka Revel'
        d['Subject'] = 'high resolution map '+rivername0
        d['Keywords'] = 'high resolution, river profile '+restag+', '+rivername0
        d['CreationDate'] = datetime.datetime(2021, 6, 24)
        d['ModDate'] = datetime.datetime.today()
    return 0
#############################
rivernames=["MEKONG"]#,"LENA","CONGO","AMAZONAS","AMAZONAS","VOLGA","MISSISSIPPI"] #"AMAZONAS"
streams=["MEKONG"]#,"LENA","CONGO","AMAZONAS","MADEIRA","VOLGA","MISSISSIPPI"] #
numrivers=len(rivernames)
#############################
inputlist=[]
for nriver in np.arange(numrivers):
    rivername0=rivernames[nriver]
    stream0=streams[nriver]
    inputlist.append([rivername0, stream0])
#==========================
#== parallel writing pdf ==
#==========================
# para_flag=1
para_flag=0
#--
if para_flag==1:
    p=Pool(6)
    list(p.map(mk_pdf,inputlist))
    p.terminate()
else:
    list(map(mk_pdf,inputlist))