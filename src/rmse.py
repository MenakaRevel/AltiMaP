#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import datetime
import numpy as np
from numpy import ma
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings;warnings.filterwarnings('ignore')
import xarray as xr
import math
import seaborn as sns
import matplotlib.patches as mpatches
import os
#========================================
#========  functions calculations  ======
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
    # return np.sqrt(np.mean((s-o)**2))
    return np.sqrt(np.ma.mean(np.ma.masked_where(o<=0.0,(s-o)**2)))
#========================================

#=========================================
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
    fname=odir+"/HydroWeb/hydroweb-hydroda_cmf_daily_wse_VIC_BC.nc"
if TAG=="CGLS":
    fname=odir+"/CGLS/cgls_cmf_daily_wse_VIC_BC.nc"
if TAG=="ICESat":
    fname=odir+"/ICESat/icesat_cmf_daily_wse_VIC_BC.nc"
if TAG=="HydroSat":
    fname=odir+"/HydroSat/hydrosat_cmf_daily_wse_VIC_BC.nc"
#==========================
fname="/cluster/data6/menaka/CaMaVal/results/HydroWeb/hydroweb_cmf_daily_wse_VIC_BC_glb_15min.nc"
nc = xr.open_dataset(fname)
sfcelv_hydroweb=nc.sfcelv_hydroweb.values
sfcelv_cmf=nc.sfcelv_cmf.values
pname=nc.name.values
lons=nc.lon.values
lats=nc.lat.values
basins=nc.Basin.values
rivers=nc.river.values
countries=nc.country.values
sfcelv_hydroweb_max=nc.sfcelv_hydroweb_max.values
sfcelv_hydroweb_min=nc.sfcelv_hydroweb_min.values
sfcelv_cmf_max=nc.sfcelv_cmf_max.values
sfcelv_cmf_min=nc.sfcelv_cmf_min.values
# disttomouth=nc.disttomouth.values
# elediff=nc.elediff.values
nc.close()
#===========================
pnum=len(pname)
#===========================
biasnum=0
rmsenum=0
#===========================
for point in np.arange(pnum):
    # if basins[point][0] != "Amazonas":
    #     continue
    rmse=RMSE(sfcelv_cmf[:,point],sfcelv_hydroweb[:,point])
    cmf_mean=np.mean(sfcelv_cmf[:,point])
    obs_mean=np.mean(ma.masked_equal(sfcelv_hydroweb[:,point],-9999.0))
    BIAS=abs(cmf_mean-obs_mean)
    # print (basins[point][0],rivers[point][0],str(pname[point][0]), " RMSE: ",rmse, "BIAS: ", BIAS)
    line="%67s%10.4f%10.4f"%(str(pname[point][0]),rmse,BIAS)
    print (line)
    # print (str(pname[point][0]),rmse,BIAS,np.mean(sfcelv_cmf[:,point]),np.mean(ma.masked_equal(sfcelv_hydroweb[:,point],-9999.0)))
    if BIAS > 5.0:
        biasnum=biasnum+1
    if rmse > 5.0:
        rmsenum=rmsenum+1
    # # print (cmf_mean,obs_mean,BIAS)
    # datayy0.append(BIAS)
    # dataxx0.append(abs(elediff[point]))

# print (float(biasnum)/float(pnum))*100.0
# print (float(rmsenum)/float(pnum))*100.0