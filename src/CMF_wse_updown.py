#!/usr/bin/env python
# coding=utf-8
# Extract CMF WSE at corresponding virtual stations and store them into a separate file
#
import os, sys
import numpy as np
from netCDF4 import * #Dataset
import datetime
import string
from numpy import ma
import xarray as xr
from scipy.interpolate import interp1d
#
# from read_discharge import read_discharge, read_discharge_multi
from read_sfcelv import read_sfcelv, read_sfcelv_multi
from read_patchMS import upstream
import read_hydroweb as hweb
import read_cgls as cgls
import read_icesat as icesat 
import read_hydrosat as hsat

curr_pwd = '/cluster/data6/menaka/CaMaVal/'

CaMa_dir = '/cluster/data6/menaka/CaMa-Flood_v396a_20200514/'

#TAG='camavali'
TAG='HydroWeb'
#TAG='CGLS'
#TAG='HydroSat'

#runoff_folder = '/cluster/data6/x.zhou/CaMa_v396/glb_06min/out/e2o_ecmwf/'
#runoff_folder = '/cluster/data6/x.zhou/CaMa_v396/glb_06min/out/e2o_ecmwf_correct/'
# runoff_folder = '/cluster/data6/x.zhou/CaMa_v396b/glb_06min/out/e2o_ecmwf_1050/'
# runoff_folder = '/cluster/data6/x.zhou/CaMa_v396/glb_06min/out/e2o_cnrs/'
#runoff_folder = '/cluster/data6/menaka/ensemble_org/CaMa_out/GLBVIC_BC001/'
runoff_folder = '/cluster/data6/menaka/ensemble_org/CaMa_out/GLBVIC_BC_USED/'

TAG=sys.argv[1]
syear=int(sys.argv[2])
eyear=int(sys.argv[3])
runoff_folder=sys.argv[4]
CaMa_dir=sys.argv[5]
curr_pwd=sys.argv[6]
out_pwd =sys.argv[7]

csize = 0.1

# out_pwd = curr_pwd+'/results/'+ TAG +'/'

# syear = 1979
# eyear = 2013

area_min=0.0 #1.0e11

peak_time_var=100

slope_threshold=100000 #m

# os.system('mkdir -p '+out_pwd)


ssize=12
#pdf = PdfPages( pwd + 'gauges.pdf')
def read_wse_multi(ix, iy, syear, eyear, number, lat, lon):
    #print ix1,iy1
    wse = np.zeros( (len(ix), nbdays), 'f')
    wse_max = np.zeros( (len(ix), nbyears), 'f')
    wse_min = np.zeros( (len(ix), nbyears), 'f')
    wse_max_loc = np.zeros( (len(ix), nbyears), 'f')
    wse_min_loc = np.zeros( (len(ix), nbyears), 'f')
    for year in range(syear, eyear+1):
        #print year
        s_days = int( (datetime.date(year , 1,1) - datetime.date(syear, 1, 1)). days)
        e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)

        f = runoff_folder + '/sfcelv'+str(year)+'.bin'
        #print f, len(ix), len(iy)
        tmp = read_sfcelv_multi( ix, iy, e_days-s_days, f, nx, ny)

        #print year, e_days - s_days, s_days, e_days, outflw.shape
        wse[:,s_days:e_days] = tmp
        wse_max[:,year-syear] = np.nanmax(tmp, axis=1)
        wse_min[:,year-syear] = np.nanmin(tmp, axis=1)
        wse_max_loc[:,year-syear] = np.argmax(tmp, axis=1)
        wse_min_loc[:,year-syear] = np.argmin(tmp, axis=1)

    return wse , wse_max, wse_min, wse_max_loc, wse_min_loc

# read WSE from CMF outputs
def read_wse(ix, iy, syear, eyear, number, lat, lon):
    wse = np.zeros( nbdays, 'f')
    wse_max = np.zeros( nbyears, 'f')
    wse_min = np.zeros( nbyears, 'f')
    wse_max_loc = np.zeros( nbyears, 'f')
    wse_min_loc = np.zeros( nbyears, 'f')
    for year in range(syear, eyear+1):
        s_days = int( (datetime.date(year , 1,1) - datetime.date(syear, 1, 1)). days)
        e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)

        f = runoff_folder + '/sfcelv'+str(year)+'.bin'

        #outflw = np.fromfile(f, 'float32').reshape(-1,len(lat_global), len(lon_global))
        #if ix2 < 0 and iy2 < 0:
        #    tmp = outflw[:,iy1, ix1]
        #else:
        #    tmp = outflw[:,iy1, ix1] + outflw[:,iy2, ix2]

        #print tmp
        #sno=len(numbers)
        tmp = read_sfcelv( ix, iy, e_days-s_days, f, nx, ny)

        #print year, e_days - s_days, s_days, e_days, outflw.shape
        wse[s_days:e_days] = tmp
        wse_max[year-syear] = np.nanmax(tmp)
        wse_min[year-syear] = np.nanmin(tmp)
        wse_max_loc[year-syear] = np.argmax(tmp)
        wse_min_loc[year-syear] = np.argmin(tmp)

        '''
        # this is to check if the location is right or not.
        fig = plt.figure(figsize=(ssize ,ssize*0.5),dpi=300)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        plt.title('GRDC: '+ str(number))

        tmp = np.nanmean(outflw, axis=0)
        tmp = np.ma.masked_where(tmp > 1.e10, tmp)
        tmp = np.ma.masked_where(tmp <=0 , tmp)
        im2=plt.imshow(tmp, cmap=cm.YlOrRd, extent=(-180,180,-90,90))

        plt.scatter( lon, lat, s= 20, c='r')
        plt.scatter( lon_global[ix1], lat_global[iy1], s=10, c='b')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.15)

        cbar=plt.colorbar(im2, cax=cax)
        cbar.set_label( 'Average discharge (m3/s)', size=ssize)
        cbar.ax.tick_params(labelsize=ssize)

        pdf.savefig()
        plt.close()

        # =====
        fig = plt.figure(figsize=(ssize ,ssize*0.5),dpi=300)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        plt.title('GRDC: '+ str(number) + ' dis(grdc):' + str(int(grdc_mean)) + ' (cmf):'+str(int(np.nanmean(dis))))

        tmp = np.nanmean(outflw, axis=0)
        tmp = np.ma.masked_where(tmp > 1.e10, tmp)
        tmp = np.ma.masked_where(tmp <=0 , tmp)
        vmax = np.nanmax( tmp[iy1-12:iy1+12, ix1-12:ix1+12] )
        im2=plt.imshow(tmp, cmap=cm.YlOrRd, extent=(-180,180,-90,90), vmax=vmax)

        plt.scatter( lon, lat, s= 20, c='r', label='GRDC')
        plt.scatter( lon_global[ix1], lat_global[iy1], s=10, c='b', label='CMF')

        ax.set_xlim( lon - 3, lon + 3)
        ax.set_ylim( lat - 3, lat + 3)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2%", pad=0.15)

        cbar=plt.colorbar(im2, cax=cax)
        cbar.set_label( 'Average discharge (m3/s)', size=ssize)
        cbar.ax.tick_params(labelsize=ssize)
        plt.legend(loc=0)

        pdf.savefig()
        plt.close()
        '''

        #print outflw.shape
        #print outflw[:,iy1,ix1]

        #print lon, lat, lon_global[ix1], lat_global[iy1]
    return wse, wse_max, wse_min, wse_max_loc, wse_min_loc
#-------------------------------------
def slope(ix,iy,nextxy,uparea,elevtn,nxtdst,rivseq):
    if rivseq[iy,ix]>1:
        nextX=nextxy[0]
        nextY=nextxy[1]
        uXX, uYY = upstream(ix+1,iy+1,nextX.T,nextY.T,uparea.T)
        uXX = uXX - 1
        uYY = uYY - 1
        slp=(elevtn[uYY,uXX]-elevtn[iy,ix])
    else:
        slp=0
    return slp
# =============================
# Read the GRDC daily data(only selected gauges)
# name = []
# river = []
# country = []
#o = Dataset('GRDC_Daily_selected.nc','r')
#o = Dataset('../shared/GRDC_Daily_selected-1958_2014.nc','r')
# o = Dataset('/cluster/data6/x.zhou/data/discharge/GRDC_2019/script/GRDC_Daily_selected-CaMaVali-1980_2014.nc','r')
# grdc_styear = 1980
# lat = o.variables['lat'][:]
# lon = o.variables['lon'][:]
# area = o.variables['area'][:]
# number = o.variables['number'][:]
# hydrographs = o.variables['hydrographs'][:]
# hydrographs[hydrographs>1.e10] = np.nan
# for i in range(len(lat)):
#     name.append(string.rstrip(string.join(o.variables['name'][i],sep='')))
#     river.append(string.rstrip(string.join(o.variables['river'][i],sep='')))
#     country.append(string.rstrip(string.join(o.variables['country'][i],sep='')))

# o.close()
#=============================
# Read the CMF variables
if csize == 0.25:
    mapname = 'glb_15min'
    nx      = 1440
    ny      = 720
elif csize == 0.1:
    mapname = 'glb_06min'
    nx      = 3600
    ny      = 1800
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
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
lonlat = np.fromfile(lonlat,np.float32).reshape(2,ny,nx)
uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)
nxtdst = np.fromfile(nxtdst,np.float32).reshape(ny,nx)
rivseq = np.fromfile(rivseq,np.int32).reshape(ny,nx)
#-----
nextX=nextxy[0]
nextY=nextxy[1]
#-----
# # if TAG=="HydroWeb":
# #     # read HydroWeb data
# #     numbers,rivers,pname,lons,lats,xlist,ylist,legm08,legm96,lsat,leled,ldist,lflag=hweb.get_hydroweb_locs(mapname)
# #     hwb_metadata=hweb.metadata()
# # elif TAG=="CGLS":
# #     # read CGLS data
# #     numbers,rivers,pname,lons,lats,xlist,ylist,legm08,legm96,lsat,leled,ldist,lflag=cgls.get_cgls_locs(mapname)
# #     hwb_metadata=cgls.metadata()
# # elif TAG=="ICESat":
# #     # read ICESat data
# #     numbers,rivers,pname,lons,lats,xlist,ylist,legm08,legm96,lsat,leled,ldist,lflag=icesat.get_icesat_locs(mapname)
# #     hwb_metadata=icesat.metadata()
# # elif TAG=="HydroSat":
# #     # read HydroSat data
# #     numbers,rivers,pname,lons,lats,xlist,ylist,legm08,legm96,lsat,leled,ldist,lflag=hsat.get_hydrosat_locs(mapname)
# #     hwb_metadata=hsat.metadata()
# # #print hwb_metadata.keys()
# # #print hwb_metadata['R_AMAZONAS_AMAZONAS_KM0344'][3]
# # # if csize == 0.25:
# # #     data = np.loadtxt('../shared/GRDC_alloc_15min.txt', skiprows=1)
# # # elif csize == 0.1:
# # #     data = np.loadtxt('../shared/GRDC_alloc_06min.txt', skiprows=1)
# # lat_global = np.arange( 90 - csize / 2., -90, -csize)
# # lon_global = np.arange( -180 + csize / 2. , 180, csize)

# number_cmf = map(int, data[:,0])
# lat_cmf = data[:,1]
# lon_cmf = data[:,2]
# # err
# area_grdc = data[:,4]
# area_cmf = data[:,5]
# diff = data[:,6]
# ups_num=map(int, data[:,7])
# ix1 = map(int, data[:,8]-1)
# iy1 = map(int, data[:,9]-1)
# ix2 = map(int, data[:,10]-1)
# iy2 = map(int, data[:,11]-1)

# area1 = data[:,12]
# area2 = data[:,13]
#==============================
# # # CaMa parameters
# # uparea = CaMa_dir+"/map/"+mapname+"/uparea.bin"
# # uparea = np.fromfile(uparea,np.float32).reshape(ny,nx)

# # subbasin = "/cluster/data6/menaka/CaMa_subbasin/output/subbasin_"+mapname+".bin"
# # subbasin = np.fromfile(subbasin,np.float32).reshape(ny,nx)
# # # =============================
# # # define the variables that will be writen to files
# # coll_area_grdc = []
# # coll_area_cmf = []

# # coll_lat = []
# # coll_lon = []
# # coll_xx = []
# # coll_yy = []
# # coll_number = []
# # coll_name = []
# # coll_river = []
# # coll_basin = []
# # coll_country = []
# # coll_subbasin = []
# # coll_satellite = []
# # coll_duration = []
# # coll_status = []
# # coll_elediff = []
# # coll_disttomouth = []
# # coll_flag = []

# # nbdays = int( (datetime.date(eyear + 1, 1,1) - datetime.date(syear, 1, 1)). days)
# # nbyears = eyear - syear + 1

# # coll_wse = np.zeros( (nbdays, len(lats), 2), 'f')
# # coll_wse_max = np.zeros( ( nbyears, len(lats), 2), 'f')
# # coll_wse_min = np.zeros( ( nbyears, len(lats), 2), 'f')



# =============================
# loop from the first gauges to check if it exists in VS definitions
# if exists, read the CMF sfcelv
count  = 0
# xxlist = []
# yylist = []
# nums   = []
# llats  = []
# llons  = []

odir=out_pwd
if TAG=="HydroWeb":
    fname0=odir+"/hydroweb_cmf_daily_wse_VIC_BC.nc"
    fname1=odir+"/hydroweb_cmf_linear_daily_wse_VIC_BC.nc"
    fname2=odir+"/hydroweb_cmf_elediff_daily_wse_VIC_BC.nc"
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

# read netCDF file
######################## original data #####################
nc0 = xr.open_dataset(fname0)
coll_number=nc0.number.values
coll_wse0=nc0.sfcelv_hydroweb.values
coll_cmf0=nc0.sfcelv_cmf.values
coll_name=nc0.name.values
coll_lon=nc0.lon.values
coll_lat=nc0.lat.values
coll_xx=nc0.variables["x-coord"].values
coll_yy=nc0.variables["y-coord"].values
coll_basin=nc0.Basin.values
coll_river=nc0.river.values
coll_country=nc0.country.values
coll_area_cmf=nc0.area_cmf.values
coll_satellite=nc0.satellite.values
coll_status=nc0.status.values
# sfcelv_hydroweb_max=nc0.sfcelv_hydroweb_max.values
# sfcelv_hydroweb_min=nc0.sfcelv_hydroweb_min.values
# sfcelv_cmf_max=nc0.sfcelv_cmf_max.values
# sfcelv_cmf_min=nc0.sfcelv_cmf_min.values
ldist=nc0.disttomouth.values
leled=nc0.elediff.values
# flag=nc0.flag.values
nc0.close()

nbdays = int( (datetime.date(eyear + 1, 1,1) - datetime.date(syear, 1, 1)). days)
nbyears = eyear - syear + 1

coll_wse = np.zeros( (nbdays, len(coll_lat), 2), 'f')
coll_wse_max = np.zeros( ( nbyears, len(coll_lat), 2), 'f')
coll_wse_min = np.zeros( ( nbyears, len(coll_lat), 2), 'f')


# get the upstream pixels
pnum=len(coll_name)
uxlist=[]
uylist=[]
for point in np.arange(pnum):
    iix=coll_xx[point]
    iiy=coll_yy[point]
    if rivseq[iiy,iix] == 1: 
        uxlist.append(iix)
        uylist.append(iiy)
    else:
        uXX, uYY = upstream(iix+1,iiy+1,nextX.T,nextY.T,uparea.T)
        uxlist.append(uXX)
        uylist.append(uYY)

dxlist=[]
dylist=[]
for point in np.arange(pnum):
    iix=coll_xx[point]
    iiy=coll_yy[point]
    if nextX[iiy,iix] == -9 or nextX[iiy,iix] == -10:
        dxlist.append(iix)
        dylist.append(iiy)
    elif nextX[iiy,iix] == -9999:
        dxlist.append(iix)
        dylist.append(iiy)
    else:
        dxlist.append(nextX[iiy,iix])
        dylist.append(nextY[iiy,iix])
# get CMF results
# wse_cmf, wse_cmf_max, wse_cmf_min, wse_cmf_max_loc, wse_cmf_min_loc = read_wse_multi(xlist, ylist, syear, eyear, numbers, llats, llons)

# get upstream of xlist and ylist
uwse_cmf, uwse_cmf_max, uwse_cmf_min, uwse_cmf_max_loc, uwse_cmf_min_loc = read_wse_multi(uxlist, uylist, syear, eyear, coll_number, coll_lat, coll_lon)

# get downstream of xlist and ylist
dwse_cmf, dwse_cmf_max, dwse_cmf_min, dwse_cmf_max_loc, dwse_cmf_min_loc = read_wse_multi(dxlist, dylist, syear, eyear, coll_number, coll_lat, coll_lon)

############
# simulations
############
# upstream
coll_wse[:, :, 0] = uwse_cmf.T

# downstream
coll_wse[:, :, 1] = dwse_cmf.T

count = len(coll_number)


# # for i in range(len(coll_number)):
# #     #print i
# #     num = coll_number[i]
    # # # print ("getting data : ",num ,pname[i])
    # # # try :
    # # # #if 1:
    # # #     a = np.where( number_cmf == number[i])[0][0]
    # # #     found = 1

    # # # except:
    # # #     found = 0
    # # #     print i, number[i]

    # # # condtion for observations
    # # #if found == 1:
    # # if TAG=="HydroWeb":
    # #     hwb_data=hweb.HydroWeb_continous_WSE(pname[i],syear,1,1,eyear,12,31,legm08[i],legm96[i])
    # # if TAG=="CGLS":
    # #     hwb_data=cgls.cgls_continous_WSE(pname[i],syear,1,1,eyear,12,31,legm08[i],legm96[i])
    # # if TAG=="ICESat":
    # #     hwb_data=icesat.ICESat_continous_WSE(pname[i],syear,1,1,eyear,12,31,legm08[i],legm96[i])
    # # if TAG=="HydroSat":
    # #     hwb_data=hsat.HydroSat_continous_WSE(pname[i],syear,1,1,eyear,12,31,legm08[i],legm96[i])
    # # #----
    # # times=len(hwb_data)
    # # observ=np.sum((hwb_data!=-9999.0)*1.0)
    # # found=np.sum((hwb_data==-9999.0)*1.0)
     
    # # #print found
    # # ################
    # # # data availabiltiy
    # # ################
    # # if found==float(times):
    # #     # print ("no data: ",pname[i])
    # #     continue
    # # ################
    # # # At least one year of data (365/35  at least around 10 observations)
    # # ################
    # # if observ<1:
    # #     # print ("no data: ",pname[i], observ)
    # #     continue

    # # ################
    # # # condtion for mainstream
    # # ################
    # # if uparea[ylist[i],xlist[i]] < area_min:
    # #     # print ("smaller river: ",pname[i], uparea[ylist[i],xlist[i]])
    # #     continue
    
    # # ################
    # # # condition for slope
    # # ################
    # # slp=slope(xlist[i],ylist[i],nextxy,uparea,elevtn,nxtdst,rivseq)
    # # if slp > slope_threshold:
    # #     # print ("high slope:", pname[i], slp)
    # #     continue

    # # # print (hwb_metadata[pname[i]][0])
    
    
    # # # for year in range(syear, eyear+1):
    # # #     s_days = int( (datetime.date(year  , 1, 1) - datetime.date(syear, 1, 1)). days)
    # # #     e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)
    # # #    coll_wse_max[year-syear,count,0] = np.nanmax( ma.masked_equal(hwb_data[s_days:e_days], -9999.0))
    # # #    coll_wse_min[year-syear,count,0] = np.nanmin( ma.masked_equal(hwb_data[s_days:e_days], -9999.0) )
    # # #     s_days_grdc = int( (datetime.date(year , 1,1) - datetime.date(grdc_styear, 1, 1)). days)
    # # #     e_days_grdc = int( (datetime.date(year+1, 1, 1) - datetime.date(grdc_styear, 1, 1)). days)

    # # #     #print s_days_grdc, e_days_grdc
    # # #     coll_dis[s_days:e_days, count, 0] = hydrographs[s_days_grdc:e_days_grdc,i]
    # # #     coll_dis_max[year-syear,count,0] = np.nanmax( hydrographs[s_days_grdc:e_days_grdc,i] )
    
    # # #coll_wse[:, count, 0] = hwb_data
    # # #hweb.HydroWeb_continous_WSE(pname[i],syear,1,1,eyear,12,31,legm08[i],legm96[i])

    # # # tmp = hydrographs[:,i]
    # # # tmp[tmp>1.e10] = np.nan

    # # # IX1.append(ix1[a])
    # # # IY1.append(iy1[a])
    # # # IX2.append(ix2[a])
    # # # IY2.append(iy2[a])
    # # #dis_cmf, dis_cmf_max = read_dis(ix1[a], iy1[a], ix2[a], iy2[a], syear, eyear, number[i], lat[i], lon[i], np.nanmean(tmp))
    # # #coll_dis[:, count, 1] = dis_cmf
    # # #coll_dis_max[:,count, 1] = dis_cmf_max

    # # #tmp1 =  np.nanmean(coll_dis[:,count,0])
    # # #tmp2 = np.nanmean(coll_dis[:,count,1])
    # # #print i, number[i], name[i], river[i], country[i], lon[i], lat[i], area[i], lon_global[ix1[a]], lat_global[iy1[a]], area_cmf[a], tmp1, tmp2, tmp2/tmp1

    # # ################
    # # # Find max min wse
    # # ################
    # # #wse_cmf, wse_cmf_max, wse_cmf_min, wse_cmf_max_loc, wse_cmf_min_loc = read_wse(xlist[i], ylist[i], syear, eyear, num, lats[i], lons[i])
    # # # # lmaxhw=[]
    # # # # lminhw=[]
    # # # # for year in np.arange(syear,eyear+1):
    # # # #     s_days = int( (datetime.date(year  , 1, 1) - datetime.date(syear, 1, 1)). days)
    # # # #     e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)
    # # # #     maxpoint=wse_cmf_max_loc[i,year-syear]
    # # # #     minpoint=wse_cmf_min_loc[i,year-syear]
    # # # #     maxrange1=int(max(s_days,s_days+maxpoint-peak_time_var))
    # # # #     maxrange2=int(min(e_days,s_days+maxpoint+peak_time_var)+1)
    # # # #     minrange1=int(max(s_days,s_days+minpoint-peak_time_var))
    # # # #     minrange2=int(max(e_days,s_days+minpoint+peak_time_var)+1)
    # # # #     maxhw=np.nanmax(ma.masked_equal(hwb_data[maxrange1:maxrange2],-9999.0).filled(-9999.0))#.filled(-9999.0)
    # # # #     minhw=np.nanmin(ma.masked_equal(hwb_data[minrange1:minrange2],-9999.0).filled(9999.0))#.filled(-9999.0)
    # # # #     #print pname[i], year ,maxpoint, minpoint, maxhw, minhw #s_days+maxpoint, s_days+minpoint
    # # # #     if maxhw != -9999.0 and minhw != 9999.0:
    # # # #         #print (point,year,maxhw,minhw)
    # # # #         lmaxhw.append(maxhw)
    # # # #         lminhw.append(minhw)
    # # # #     else:
    # # # #         lmaxhw.append(-9999.0)
    # # # #         lminhw.append(-9999.0)

    # # # # ###########
    # # # # # condition for high-low flow
    # # # # ###########
    # # # # if ( np.sum((np.array(lmaxhw)!=-9999.0)*1.0) == 0.0 or np.sum((np.array(lminhw)!=-9999.0)*1.0) == 0.0 ):
    # # # #     # print ("no high-low flow")
    # # # #     continue

    # # # print ("Included : ",num ,pname[i])

    # # ###########
    # # # append data
    # # ###########
    # # xxlist.append(xlist[i])
    # # yylist.append(ylist[i])
    # # nums.append(num)
    # # llats.append(lats[i])
    # # llons.append(lons[i])
    # # coll_xx.append(xlist[i])
    # # coll_yy.append(ylist[i])
    # # # coll_area_grdc.append(uparea[ylist[i],xlist[i]])
    # # coll_area_cmf.append(uparea[ylist[i],xlist[i]])
    # # coll_lat.append(lats[i])
    # # coll_lon.append(lons[i])
    # # coll_number.append(numbers[i])
    # # coll_name.append(pname[i])
    # # #print hwb_metadata[pname[i]]
    # # #print hwb_metadata[pname[i]][0]
    # # coll_river.append(hwb_metadata[pname[i]][0])
    # # coll_basin.append(hwb_metadata[pname[i]][1])
    # # coll_country.append(hwb_metadata[pname[i]][2])
    # # coll_subbasin.append(subbasin[ylist[i],xlist[i]])
    # # coll_satellite.append(hwb_metadata[pname[i]][3])
    # # # coll_duration.append(hwb_metadata[pname[i]][4]+"-"+hwb_metadata[pname[i]][5])
    # # coll_status.append(hwb_metadata[pname[i]][6])
    # # # coll_elediff.append(leled[i])
    # # # coll_disttomouth.append(ldist[i])
    # # # coll_flag.append(lflag[i])

    # # ########################
    # # # observed high low flows
    # # ########################
    # # # coll_wse_max[:,count,0]=np.array(lmaxhw)
    # # # coll_wse_min[:,count,0]=np.array(lminhw)

    # # ############
    # # # observations
    # # ############
    # # coll_wse[:, count, 0] = hwb_data

# #     ############
# #     # CMF data
# #     ############
# #     if rivseq[coll_yy[i],coll_xx[i]] == 1:
# #         total_dist=1.0
# #         coll_wse[:,count,1] = coll_cmf0[:,i]
# #     else:
# #         total_dist=nxtdst[uylist[i],uxlist[i]]*1e-3
# #         # lldist=np.array([0.0,total_dist])
# #         # llcmf=np.array([uwse_cmf.T[:,i],wse_cmf.T[:,i]])
# #         # interfunc=interp1d(lldist,llcmf)
# #         coll_wse[:,count,1] = coll_cmf0[:,i] + ((uwse_cmf.T[:,i]-coll_cmf0[:,i])/total_dist)*ldist[i]
# #         # coll_wse[:,count,1] = interfunc(ldist[i])
# #     # coll_wse_max[:,count,1] = wse_cmf_max.T[:,i]
# #     # coll_wse_min[:,count,1] = wse_cmf_min.T[:,i] 
    
# #     print (coll_name[i][0], rivseq[coll_yy[i],coll_xx[i]], nxtdst[uylist[i],uxlist[i]]*1e-3, ldist[i], np.mean(coll_wse[:,count,1]))       
# #     count += 1

# # print (count)
# wse_cmf, wse_cmf_max, wse_cmf_min, wse_cmf_max_loc, wse_cmf_min_loc = read_wse_multi(xxlist, yylist, syear, eyear, nums, llats, llons)
# coll_wse[:,:count,1] = wse_cmf.T
# coll_wse_max[:,:count,1] = wse_cmf_max.T
# coll_wse_min[:,:count,1] = wse_cmf_min.T

# ########
# for point in np.arange(count):
#     for year in np.arange(syear,eyear+1):
#         s_days = int( (datetime.date(year  , 1, 1) - datetime.date(syear, 1, 1)). days)
#         e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)
#         maxpoint=wse_cmf_max_loc[point,year-syear]
#         minpoint=wse_cmf_min_loc[point,year-syear]
#         maxrange1=int(max(s_days,s_days+maxpoint-peak_time_var))
#         maxrange2=int(min(e_days,s_days+maxpoint+peak_time_var)+1)
#         minrange1=int(max(s_days,s_days+minpoint-peak_time_var))
#         minrange2=int(max(e_days,s_days+minpoint+peak_time_var)+1)
#         maxhw=np.nanmax(ma.masked_equal(coll_wse[maxrange1:maxrange2,point,0],-9999.0).filled(-9999.0))#.filled(-9999.0)
#         minhw=np.nanmin(ma.masked_equal(coll_wse[minrange1:minrange2,point,0],-9999.0).filled(9999.0))#.filled(-9999.0)
#         #print pname[point][0], year ,maxsfcelv, minsfcelv, maxhw, minhw #s_days+maxpoint, s_days+minpoint
#         if maxhw != -9999.0 and minhw != 9999.0:
#             #print (point,year,maxhw,minhw)
#             coll_wse_max[syear-year,point,0]=maxhw
#             coll_wse_min[syear-year,point,0]=minhw
#         else:
#             coll_wse_max[syear-year,point,0]=-9999.0
#             coll_wse_min[syear-year,point,0]=-9999.0
# write everything to the file
#o = Dataset(out_pwd + 'grdc_cmf_daily_dis_multi.nc', 'w')
#o = Dataset(out_pwd + 'grdc_cmf_daily_dis_correct.nc', 'w')
#o = Dataset(out_pwd + 'grdc_cmf_daily_dis_correct_roughness.nc', 'w')
#o = Dataset(out_pwd + 'grdc_cmf_daily_dis_multi_cnrs.nc', 'w')

if TAG=="HydroWeb":
    fname='hydroweb_cmf_updown_daily_wse_VIC_BC.nc'
if TAG=="CGLS":
    fname='cgls_cmf_updown_daily_wse_VIC_BC.nc'
if TAG=="ICESat":
    fname='icesat_cmf_updown_daily_wse_VIC_BC.nc'
if TAG=="HydroSat":
    fname='hydrosat_cmf_updown_daily_wse_VIC_BC.nc'
#-------
print (out_pwd+"/"+fname, count)
o = Dataset(out_pwd+"/"+fname, 'w')
#o = Dataset(out_pwd + 'cgls_cmf_daily_wse_VIC_BC.nc', 'w')
o.history='Daily water surface elevation starts from '+str(syear)+'-01-01, The water surface elevations are in WGS84 EGM96\ncreated by Menaka Revel @IIS, U-Tokyo'
#o.threshold='Larger than 10000km2 and with >15 years data in '+str(syear)+'-2014'
o.threshold='data in '+str(syear)+'-'+str(eyear)

o.createDimension('lat', len(coll_lat))
o.createDimension('lon', len(coll_lon))
o.createDimension('stations', len(coll_lat))
o.createDimension('nbdays', nbdays)
o.createDimension('nbyears', nbyears)
strlen=1
o.createDimension('str', strlen)

v=o.createVariable('lat', 'f', ('lat',))
v.longname='latitude'
v[:] = coll_lat[:]

v=o.createVariable('lon', 'f', ('lon',))
v.longname='longitude'
v[:] = coll_lon[:]

# v=o.createVariable('y-coord', 'i', ('lat',))
# v.longname='y coordinate of the CMF grid'
# v[:] = coll_yy[:]

# v=o.createVariable('x-coord', 'i', ('lon',))
# v.longname='x coordinate of the CMF grid'
# v[:] = coll_xx[:]

stname=o.createVariable('name', 'S60', ('stations', 'str',))
stname.long_name="Virtual Station Name"
stname[:] = coll_name 
# str_out   = stringtochar(np.array(coll_name)) #, dtype='object')
# stname[:] = str_out
# for i in range(count) :
#     nana=coll_name[i]
#     stname[i,:]=chartostring(np.asarray(list(nana.ljust(strlen))))

riv=o.createVariable('river', 'S60', ('stations','str',))
riv.long_name="River Name"
riv[:] = coll_river
# str_out= stringtochar(np.array(coll_river))  #, dtype='object')
# riv[:] = str_out
# for i in range(count) :
#     nana=coll_river[i]
#     riv[i,:]=chartostring(np.asarray(list(nana.ljust(strlen))))

bsn=o.createVariable('Basin', 'S60', ('stations','str',))
bsn.long_name="Basin Name"
bsn[:] = coll_basin
# str_out= stringtochar(np.array(coll_river))  #, dtype='object')
# bsn[:] = str_out
# for i in range(count) :
#     nana=coll_basin[i]
#     bsn[i,:]=chartostring(np.asarray(list(nana.ljust(strlen))))


ctry=o.createVariable('country', 'S60', ('stations','str',))
ctry.long_name="Country of station"
ctry[:] = coll_country
# str_out = stringtochar(np.array(coll_country)) #, dtype='object'))
# ctry[:] = str_out
# for i in range(count) :
#     nana=coll_country[i]
#     ctry[i,:]=chartostring(np.asarray(list(nana.ljust(strlen))))

# dura=o.createVariable('duration', 'S60', ('stations','str',))
# dura.long_name="Duration of data availability"
# # str_out = stringtochar(np.array(coll_duration)) #, dtype='object'))
# # dura[:] = str_out
# for i in range(count) :
#     nana=coll_duration[i]
#     dura[i,:]=chartostring(np.asarray(list(nana.ljust(strlen))))

sat=o.createVariable('satellite', 'S60', ('stations','str',))
sat.long_name="Names of the satellites"
sat [:] = coll_satellite
# str_out = stringtochar(np.array(coll_duration)) #, dtype='object'))
# sat [:] = str_out
# for i in range(count) :
#     nana=coll_satellite[i]
#     sat[i,:]=chartostring(np.asarray(list(nana.ljust(strlen))))

stat=o.createVariable('status', 'S60', ('stations','str',))
stat.long_name="Status of the satellite"
stat [:] = coll_status
# str_out = stringtochar(np.array(coll_duration)) #, dtype='object'))
# stat [:] = str_out
# for i in range(count) :
#     nana=coll_status[i]
#     stat[i,:]=chartostring(np.asarray(list(nana.ljust(strlen))))

# v=o.createVariable('area_vs', 'f', ('stations',))
# v.longname='upstream area of virtual station'
# v.unit='km2'
# v[:] = coll_area_grdc[:]

v=o.createVariable('area_cmf', 'f', ('stations',))
v.longname='upstream area of CMF'
v.unit='km2'
v[:] = coll_area_cmf[:]

v=o.createVariable('number', 'i', ('stations',))
v.longname='Vrtual Station ID'
v[:] = coll_number[:]

# v=o.createVariable('subbasin', 'f', ('stations',))
# v.longname='CMF sub basin ID'
# v[:] = coll_subbasin[:]

# v=o.createVariable('elediff', 'f', ('stations',))
# v.longname='elevation differnce between VS and unit-catchment outlet [m]'
# v[:] = coll_elediff[:]

# v=o.createVariable('disttomouth', 'f', ('stations',))
# v.longname='distance to unit-catchment outlet from VS [km]'
# v[:] = coll_disttomouth[:]

# v=o.createVariable('flag', 'i', ('stations',))
# v.longname='Flag for VS allocation | 1-river channel | 2-unit-catchment outlet | 3-correction form land grid | 4-correction from ocean grid | 5-braided rivers'
# v[:] = coll_flag[:]

v=o.createVariable('sfcelv_upstream_cmf', 'f', ('nbdays', 'stations'))
v.longname='water surface elevation of the CMF estimation at upstream pixel'
for i in range(count):
    v[:,i] = coll_wse[:,i,0]

v=o.createVariable('sfcelv_downstream_cmf', 'f', ('nbdays', 'stations'))
v.longname='water surface elevation of the CMF estimation at downstream pixel'
for i in range(count):
    v[:,i] = coll_wse[:,i,1]

# v=o.createVariable('sfcelv_hydroweb_max', 'f', ('nbyears', 'stations'))
# v.longname='Annual maximum water surface elevation at the Vrtual Stations'
# for i in range(count):
#     v[:,i] = coll_wse_max[:,i,0]

# v=o.createVariable('sfcelv_hydroweb_min', 'f', ('nbyears', 'stations'))
# v.longname='Annual minimum water surface elevation at the Vrtual Stations'
# for i in range(count):
#     v[:,i] = coll_wse_min[:,i,0]

# v=o.createVariable('sfcelv_cmf_max', 'f', ('nbyears', 'stations'))
# v.longname='Annual maximum water surface elevation at the CMF estimation'
# for i in range(count):
#     v[:,i] = coll_wse_max[:,i,1]

# v=o.createVariable('sfcelv_cmf_min', 'f', ('nbyears', 'stations'))
# v.longname='Annual minimum water surface elevation at the CMF estimation'
# for i in range(count):
#     v[:,i] = coll_wse_min[:,i,1]

o.close()

#pdf.close()