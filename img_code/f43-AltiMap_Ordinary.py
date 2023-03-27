#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm,Normalize,ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import sys
import os
import errno
from numpy import ma 
import matplotlib.gridspec as gridspec
import string
from scipy.fftpack import fft, ifft, fftfreq
#from slacker import Slacker
from multiprocessing import Pool
from multiprocessing import Process
import xarray as xr
from scipy.stats.kde import gaussian_kde
from matplotlib.mlab import find
import pandas as pd 
import seaborn as sns
import math
import my_colorbar as mbar
import cartopy.crs as ccrs
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import re
import xarray as xr
from matplotlib.mlab import find
import seaborn as sns
#
# from read_patchMS import upstream
#==============================
CaMa_dir = "/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
mapname= "glb_06min"
#==============================
#====
# sfcelv
# odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
odir1 = '/cluster/data6/menaka/Altimetry/results/HydroWeb'
fname1 = odir1+"/hydroweb_cmf_daily_wse_VIC_BC.nc"
nc1 = xr.open_dataset(fname1)
sfcelv_hydroweb1=nc1.sfcelv_hydroweb.values
sfcelv_cmf1=nc1.sfcelv_cmf.values
lons1=nc1.lon.values
lats1=nc1.lat.values
pnames1=nc1.name.values
pnum1=len(pnames1)
print np.shape(sfcelv_cmf1), pnum1
#- masked out -9999 values
sfcelv_cmf1=ma.masked_where(sfcelv_hydroweb1==-9999.0,sfcelv_cmf1).filled(-9999.0)
#sfcelv_hydroweb=ma.masked_where(sfcelv_hydroweb==-9999.0,sfcelv_hydroweb)
sfcelv_diff1=ma.masked_where(sfcelv_hydroweb1==-9999.0,(sfcelv_cmf1-sfcelv_hydroweb1)**2).filled(-9999.0)
sfcelv_rmse1=np.mean(ma.masked_less_equal(sfcelv_diff1,0.0),axis=0)#.compressed()#
sfcelv_rmse1=sfcelv_rmse1.filled()
print np.shape(sfcelv_rmse1), type(sfcelv_rmse1)
print sfcelv_rmse1#[0:10]
# sfcelv_bias_com=ma.masked_equal(sfcelv_bias,-9999.0).compressed()


# ====
# sfcelv
# odir = "/cluster/data6/menaka/CaMaVal/results_daily/camavali"
odir2 = '/cluster/data6/menaka/Altimetry/results/HydroWeb'
fname2 = odir2+"/hydroweb_cmf_daily_wse_VIC_BC_ordinary.nc"
nc2 = xr.open_dataset(fname2)
sfcelv_hydroweb2=nc2.sfcelv_hydroweb.values
sfcelv_cmf2=nc2.sfcelv_cmf.values
lons2=nc2.lon.values
lats2=nc2.lat.values
pnames2=nc2.name.values
pnum2=len(pnames2)
print np.shape(sfcelv_cmf2), pnum2
#- masked out -9999 values
sfcelv_cmf2=ma.masked_where(sfcelv_hydroweb2==-9999.0,sfcelv_cmf2).filled(-9999.0)
#sfcelv_hydroweb=ma.masked_where(sfcelv_hydroweb==-9999.0,sfcelv_hydroweb)
sfcelv_diff2=ma.masked_where(sfcelv_hydroweb2==-9999.0,(sfcelv_cmf2-sfcelv_hydroweb2)**2).filled(-9999.0)
sfcelv_rmse2=np.mean(ma.masked_equal(sfcelv_diff2,-9999.0),axis=0)#.compressed()#
sfcelv_rmse2=sfcelv_rmse2.filled()
print np.shape(sfcelv_rmse2)
print sfcelv_rmse2[0:10]
# sfcelv_bias_com=ma.masked_equal(sfcelv_bias1,-9999.0).compressed()
