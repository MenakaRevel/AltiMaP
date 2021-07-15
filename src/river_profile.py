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
    return AIC
#==================
mapname=sys.argv[1] #"glb_06min"
restag=sys.argv[2] #"3sec"
noobs=sys.argv[3]
obstxt=sys.argv[4] #"/cluster/data6/menaka/Altimetry/fig/river_profile"
CaMa_dir=sys.argv[5] #"/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
ele_thr=0.35 # sys.argv[6] # 1.0 m observation error can also be considered. 
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
# # # rivername0="CONGO" #"AMAZONAS" #"AMAZONAS"
# # # stream0="CONGO" #"AMAZONAS" #
# # # dataname="HydroWeb"
# # odir="/cluster/data6/menaka/Altimetry/fig/river_profile"
# # # TAG="HydroWeb"
# # mapname="glb_06min"
# # restag="3sec"
# # CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
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
#======================================================================
noVS=[]
with open(noobs,"r") as f:
    lines=f.readlines()
    for line in lines:
        line    = filter(None,re.split(" ",line))
        station = line[0].strip()
        noVS.append(station)
#======================================================================
# rivername0="CONGO" #"AMAZONAS" #"AMAZONAS"
# stream0="CONGO"
fname=obstxt
############################################################
with open(fname,"r") as f:
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        #print line
        num     = line[0]
        station = line[1].strip()
        line2   = re.split("_",station)
        riv     = line2[1]
        stream  = line2[2]
        lon     = float(line[3])
        lat     = float(line[4])
        ix      = int(line[5])-1
        iy      = int(line[6])-1
        elem    = float(line[7])
        EGM08   = float(line[8])
        EGM96   = float(line[9])
        sat     = line[10].strip()
        dist    = float(line[11])
        flag    = int(line[12])
        kx      = int(line[13])
        ky      = int(line[14])
        #=======================
        # if riv != rivername0:
        #     continue
        # if riv == "AMAZONAS":
        #     if stream != stream0:
        #         continue
        # elif riv == "CONGO":
        #     if stream != stream0:
        #         continue
        # else:
        #     if stream != stream0:
        #         continue
        #-----------------------
        if station in noVS:
            continue
        #=====================================
        # print (station)
        if abs(elem - elevtn[iy,ix]) < ele_thr:
            print "%69s%5d%5d%7.3f%7.3f%7.3f"%(station,0,0,elem - elevtn[iy,ix],0.0,0.0)
        else:   #elev,meanW,meanWSE_VICBC[iYY,iXX])
            # print (lon, lat)
            west, south = westsouth(lat,lon)
            # print (west, south)
            north = south + 10.0
            east  = west + 10.0
            cname0 = cname(lat,lon)
            #========================================================
            # high-resolution data
            # print (cname0)
            visual=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".visual.bin"
            # print (visual)
            visual=np.fromfile(visual,np.int8).reshape(12000,12000)
            elev1m=CaMa_dir+"/map/"+mapname+"/"+restag+"/"+cname0+".elevtn.bin"
            # print (elev1m)
            elev1m=np.fromfile(elev1m,np.float32).reshape(12000,12000)
            #---------------------------------------------------------
            try:
                # print ("calculate river profile...........",kx,ky) #,visual[ky-1,kx-1])
                length, elevation = river_along_profile(kx,ky,west,south,res,nx,ny,hiresmap)
                #====================================
                xf = np.array([length[0],length[-1]])
                yf = np.array([elevation[0],elevation[-1]])
                M  = 1e5
                lrsq=[]
                lAIC=[]
                # print ("Fitting the polynomial..........")
                for n in np.arange(1,5+1):
                    try:
                        params = cpf(length, elevation, xf, yf, n, M)
                        # print (params)
                        r_squared = rsq(elevation,params(length))
                        AICval = AIC_val(elevation,params(length),n)
                        lrsq.append(r_squared)
                        lAIC.append(AICval)
                        # print (r_squared, AICval)
                        # rsq_char="$r^2$ : %3.2f"%(r_squared)
                        # AIC_char="$A\^IC$ : %3.2f"%(AICval)
                        # print (n, r_squared, AICval)
                    except:
                        lrsq.append(-9999)
                        lAIC.append(1e20)
                lrsq=np.array(lrsq)
                lAIC=np.array(lAIC)
                # replace nan values
                lrsq=np.nan_to_num(lrsq, nan=-9999)
                lAIC=np.nan_to_num(lAIC, nan=1e20)
                #------------------
                AIC_min=np.min(lAIC)
                AIC_mloc=np.argmin(lAIC)
                # compare with linear treand
                deltaAIC=lAIC[0]-AIC_min
                if deltaAIC < 2:
                    if lrsq[0] > 0.8:
                        # print station, deltaAIC, lrsq[0]
                        print "%69s%5d%5d%7.3f%7.3f%7.3f"%(station,1,1,elem - elevtn[iy,ix],lrsq[0],lAIC[0])
                    else:
                        print "%69s%5d%5d%7.3f%7.3f%7.3f"%(station,2,AIC_mloc+1,elem - elevtn[iy,ix],lrsq[AIC_mloc],lAIC[AIC_mloc])
                else:
                    print "%69s%5d%5d%7.3f%7.3f%7.3f"%(station,2,AIC_mloc+1,elem - elevtn[iy,ix],lrsq[AIC_mloc],lAIC[AIC_mloc])
            except:
                print "%69s%5d%5d%7.3f%7.3f%7.3f"%(station,-9,0,elem - elevtn[iy,ix],0.0,0.0)