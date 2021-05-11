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
odir="/cluster/data6/menaka/Altimetry/results"
if TAG=="HydroWeb":
    fname=odir+"/HydroWeb/hydroweb_cmf_elediff_daily_wse_VIC_BC.nc"
if TAG=="CGLS":
    fname=odir+"/CGLS/cgls_cmf_elediff_daily_wse_VIC_BC.nc"
if TAG=="ICESat":
    fname=odir+"/ICESat/icesat_elediff_cmf_daily_wse_VIC_BC.nc"
if TAG=="HydroSat":
    fname=odir+"/HydroSat/hydrosat_cmf_elediff_daily_wse_VIC_BC.nc"
nc = xr.open_dataset(fname)
sfcelv_hydroweb=nc.sfcelv_hydroweb.values
sfcelv_cmf=nc.sfcelv_cmf.values
pname=nc.name.values
lons=nc.lon.values
lats=nc.lat.values
basins=nc.Basin.values
rivers=nc.river.values
countries=nc.country.values
# sfcelv_hydroweb_max=nc.sfcelv_hydroweb_max.values
# sfcelv_hydroweb_min=nc.sfcelv_hydroweb_min.values
# sfcelv_cmf_max=nc.sfcelv_cmf_max.values
# sfcelv_cmf_min=nc.sfcelv_cmf_min.values
# disttomouth=nc.disttomouth.values
# elediff=nc.elediff.values
# flag=nc.flag.values
nc.close()
pnum=len(pname)
#print np.shape(sfcelv_hydroweb)
colors=['xkcd:pastel blue','xkcd:aqua green']
labels=["CaMa-Flood",TAG]
#==========
#nbyears = eyear - syear + 1
#nt, pnum = np.shape(sfcelv_cmf)
# print pnum#, nt
# max_hwb=np.zeros( (nbyears, pnum), 'f')
# min_hwb=np.zeros( (nbyears, pnum), 'f')
# max_cmf=np.zeros( (nbyears, pnum), 'f')
# min_cmf=np.zeros( (nbyears, pnum), 'f')
# # print (nc) 
# # print nc.group_by("time.year").mean()
# for point in np.arange(pnum):
#     for year in np.arange(syear,eyear+1):
#         s_days = int( (datetime.date(year  , 1, 1) - datetime.date(syear, 1, 1)). days)
#         e_days = int( (datetime.date(year+1, 1, 1) - datetime.date(syear, 1, 1)). days)
#         #rint year, start_dt, end_dt
#         #--------
#         maxsfcelv=np.amax(sfcelv_cmf[s_days:e_days,point])
#         minsfcelv=np.amin(sfcelv_cmf[s_days:e_days,point])
#         maxpoint=np.argmax(sfcelv_cmf[s_days:e_days,point])
#         minpoint=np.argmin(sfcelv_cmf[s_days:e_days,point])
#         maxrange1=max(s_days,s_days+maxpoint-15)
#         maxrange2=min(e_days,s_days+maxpoint+15)+1
#         minrange1=max(s_days,s_days+minpoint-15)
#         minrange2=max(s_days,s_days+minpoint+15)+1
#         maxhw=np.nanmax(ma.masked_equal(sfcelv_hydroweb[maxrange1:maxrange2,point],-9999.0).filled(-9999))#.filled(-9999.0)
#         minhw=np.nanmin(ma.masked_equal(sfcelv_hydroweb[minrange1:minrange2,point],-9999.0).filled(9999))#.filled(-9999.0)
#         #print pname[point][0], year ,maxsfcelv, minsfcelv, maxhw, minhw #s_days+maxpoint, s_days+minpoint
#         if maxhw != -9999.0 and minhw != 9999.0:
#             max_hwb[syear-year, point] = maxhw
#             min_hwb[syear-year, point] = minhw
#             max_cmf[syear-year, point] = maxsfcelv
#             min_cmf[syear-year, point] = minsfcelv
#         else:
#             max_hwb[syear-year, point] = -9999.0
#             min_hwb[syear-year, point] = -9999.0
#             max_cmf[syear-year, point] = -9999.0
#             min_cmf[syear-year, point] = -9999.0
# #----------
# max_hwb = np.array(max_hwb)
# min_hwb = np.array(min_hwb)
# max_cmf = np.array(max_cmf)
# min_cmf = np.array(min_cmf)
#==========
if TAG=="HydroWeb":
    pdfname=odir+"/HydroWeb/hydroweb_cmf_elediff_daily_wse_VIC_BC.pdf"
if TAG=="CGLS":
    pdfname=odir+"/CGLS/cgls_cmf_elediff_daily_wse_VIC_BC.pdf"
if TAG=="ICESat":
    pdfname=odir+"/ICESat/icesat_cmf_elediff_daily_wse_VIC_BC.pdf"
if TAG=="HydroSat":
    pdfname=odir+"/HydroSat/hydrosat_cmf_elediff_daily_wse_VIC_BC.pdf"
#============================
with PdfPages(pdfname) as pdf:
    for point in np.arange(pnum):
        #prepare
        cmf = sfcelv_cmf[:,point]
        org=ma.masked_where(sfcelv_hydroweb[:,point]==-9999.0,sfcelv_hydroweb[:,point])
        locs=np.where(sfcelv_hydroweb[:,point]!=-9999.0)[0]
        org=org.compressed()
        hgt=11.69
        wdt=8.27
        fig=plt.figure(figsize=(wdt, hgt))
        #plt.title(pname[point][0],fontsize=12)
        G = gridspec.GridSpec(3,2)
        ax0 = fig.add_subplot(G[0,:])
        #ax0.set_title(pname[point][0],fontsize=14)
        ax0.text(0.5,1.3,pname[point][0],va="center",ha="center",transform=ax0.transAxes,fontsize=14)
        print (point, pname[point][0])
        chctry="Country: %s"%(countries[point][0])
        ax0.text(0.1,1.1,chctry,transform=ax0.transAxes,fontsize=12)
        chbsn="Basin: %s"%(basins[point][0])
        ax0.text(0.1,1.0,chbsn,transform=ax0.transAxes,fontsize=12)
        chriv="River: %s"%(basins[point][0])
        ax0.text(0.1,0.9,chriv,transform=ax0.transAxes,fontsize=12)
        chlon="lon : %5.2f"%(lons[point])
        chlat="lat : %5.2f"%(lats[point])
        ax0.text(0.1,0.8,chlon,transform=ax0.transAxes,fontsize=12)
        ax0.text(0.1,0.7,chlat,transform=ax0.transAxes,fontsize=12)
        # 2nd column 
        # ele_dff="elevation differnce : %5.2f"%(elediff[point])
        # disttom="distance to mouth : %5.2f"%(disttomouth[point])
        # ax0.text(0.1,0.6,ele_dff,transform=ax0.transAxes,fontsize=12)
        # ax0.text(0.1,0.5,disttom,transform=ax0.transAxes,fontsize=12)
        # meancmf="mean (CaMa-Flood) : %5.2f"%(np.mean(cmf))
        # meanhwb="mean ("+TAG+")   : %5.2f"%(np.mean(org))
        # ax0.text(0.1,0.6,meancmf,transform=ax0.transAxes,fontsize=12)
        # ax0.text(0.1,0.5,meanhwb,transform=ax0.transAxes,fontsize=12)
        meanbias="mean bias : %5.2f"%(np.mean(cmf)-np.mean(org))
        ax0.text(0.1,0.6,meanbias,transform=ax0.transAxes,fontsize=12)
        RMSE1=RMSE(sfcelv_cmf[:,point],sfcelv_hydroweb[:,point])
        meanbias="RMSE : %5.2f"%(RMSE1)
        ax0.text(0.1,0.5,meanbias,transform=ax0.transAxes,fontsize=12)
        # amplitude
        # HydroWeb
        # mean_max_hwb=np.ma.mean(ma.masked_equal(sfcelv_hydroweb_max[:,point],-9999.0))#.filled(-9999.0)
        # mean_min_hwb=np.mean(ma.masked_equal(sfcelv_hydroweb_min[:,point],-9999.0))
        # if (np.sum((sfcelv_hydroweb_max[:,point]!=-9999.0)*1.0)==0.0 or np.sum((sfcelv_hydroweb_min[:,point]!=-9999.0)*1.0)==0.0):
        #     print ("**** no max or min ****", pname[point][0])
        #     mean_max_hwb=0.0
        #     mean_min_hwb=0.0
        # #print mean_max_hwb, np.isnan(mean_max_hwb), mean_min_hwb, np.isnan(mean_min_hwb)
        # ch_mean_max_hwb="maximum ("+TAG+") : %5.2f"%(np.nan_to_num(mean_max_hwb))
        # ch_mean_min_hwb="minimum ("+TAG+") : %5.2f"%(np.nan_to_num(mean_min_hwb))
        # ax0.text(0.1,0.3,ch_mean_max_hwb,transform=ax0.transAxes,fontsize=12)
        # ax0.text(0.1,0.2,ch_mean_min_hwb,transform=ax0.transAxes,fontsize=12)
        # CaMa-Flood
        # # mean_max_cmf=np.mean(ma.masked_equal(sfcelv_cmf_max[:,point],-9999.0))
        # # mean_min_cmf=np.mean(ma.masked_equal(sfcelv_cmf_min[:,point],-9999.0))
        # # # ch_mean_max_cmf="maximum (CaMa-Flood) : %5.2f"%(np.nan_to_num(mean_max_cmf))
        # # # ch_mean_min_cmf="minimum (CaMa-Flood) : %5.2f"%(np.nan_to_num(mean_min_cmf))
        # # # ax0.text(0.1,0.1,ch_mean_max_cmf,transform=ax0.transAxes,fontsize=12)
        # # # ax0.text(0.1,0.0,ch_mean_min_cmf,transform=ax0.transAxes,fontsize=12)
        # # amp_mean_diff=(mean_max_hwb-mean_min_hwb)-(mean_max_cmf-mean_min_cmf)
        # # #ch_amx_hwb=np.mean(ma.masked_where(max_hwb[:,point]==-9999.0,max_hwb[:,point]).compressed())
        # # meanamp="amplitiude differnce : %5.2f"%(amp_mean_diff)
        # # ax0.text(0.1,0.3,meanamp,transform=ax0.transAxes,fontsize=12)
        # # # ax0.text(0.1,-0.1,meanamp,transform=ax0.transAxes,fontsize=12)
        # # flagch="flag : %d"%(flag[point])
        # # ax0.text(0.1,0.2,flagch,transform=ax0.transAxes,fontsize=12)
        #----------------
        ax0.set_axis_off()
        ax0.spines['top'].set_visible(False)
        ax0.spines['right'].set_visible(False)
        ax0.spines['left'].set_visible(False)
        ax0.spines['bottom'].set_visible(False)
        #ax0.outline_patch.set_linewidth(0.0)
        ###########
        ax1 = fig.add_subplot(G[1,:])
        #cmf
        lines=[ax1.plot(np.arange(0,len(sfcelv_hydroweb[:,point])),cmf,color="xkcd:cornflower",label="CaMa-Flood",linewidth=0.5)[0]]
        #hydroweb
        lines.append(ax1.plot(locs,org,color="k",label=TAG,linestyle='None',linewidth=0,marker="o",fillstyle="none",markersize=5)[0])
        ax1.set_ylabel('WSE $(m)$', color='k',fontsize=10)
        ax1.tick_params('y',labelsize=8, colors='k')
        ax1.set_xlabel('Years', color='k',fontsize=10)
        ax1.tick_params('x',labelsize=8, colors='k')
        #ax1.set_title(pname[point][0],fontsize=12)
        #ax1.text(0.01,0.90,"WSE "+pname[point][0],transform=ax1.transAxes,fontsize=8)
        ax1.set_xlim(xmin=0,xmax=len(sfcelv_hydroweb[:,point])+1)
        xxlab=np.arange(syear,eyear+1,5)
        dt=int(math.ceil(((eyear-syear)+2)/5.0))
        xxlist=np.linspace(0,len(sfcelv_hydroweb[:,point]),dt,endpoint=True)
        ax1.set_xticks(xxlist)
        ax1.set_xticklabels(xxlab,fontsize=8)
        plt.legend(lines,labels,ncol=2,loc='upper right',bbox_to_anchor=(1.0, 1.13))
        ##########
        #boxplot
        ax2 = fig.add_subplot(G[2,0])
        flierprops = dict(marker='o', markerfacecolor='none', markersize=12,linestyle='none', markeredgecolor='k')
        boxprops = dict(color='grey')#facecolor='none'
        whiskerprops = dict(color='grey',linestyle="--")
        capprops = dict(color='grey')
        medianprops = dict(color='r')
        box=ax2.boxplot([cmf,org],labels=labels,boxprops=boxprops,showfliers=False, 
                        whiskerprops=whiskerprops,capprops=capprops,medianprops=medianprops, 
                        notch=False, sym=None, vert=True, whis=1.5,positions=None, widths=None, 
                        patch_artist=True,bootstrap=None, usermedians=None, conf_intervals=None)#flierprops=flierprops,
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)
        ax2.set_ylabel('WSE $(m)$', color='k',fontsize=10)
        ax2.tick_params('y',labelsize=8, colors='k')
        ax2.tick_params('x',labelsize=8, colors='k')#,labelrotation=45)
        ax2.set_xticklabels(labels,rotation=0)
        ##########
        #pdf
        ax3 = fig.add_subplot(G[2,1])
        sns.distplot(cmf, ax=ax3, hist=True, color=colors[0], label="CaMa-Flood") #ax=ax3,
        sns.distplot(org, ax=ax3, hist=True, color=colors[1], label=TAG) #ax=ax3,
        ax3.set_ylabel('density', color='k',fontsize=10)
        ax3.tick_params('y',labelsize=6, colors='k')
        ax3.set_xlabel('WSE $(m)$', color='k',fontsize=10)
        ax3.tick_params('x',labelsize=6, colors='k')
        #ax3.set_title("Histogram of Bias",fontsize=8)
        #ax3.set_xlim(xmin=-20.0,xmax=20.0)
        #ax3.text(0.01,0.95,"b",transform=ax3.transAxes,fontsize=8)
        ax3.legend(ncol=2,bbox_to_anchor=(1.0, -0.12), loc=1)
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()
    # set the file's metadata via the PdfPages object:
    d = pdf.infodict()
    d['Title'] = TAG+' and VIC BC comparison [elevation difference'
    d['Author'] = 'Menaka Revel'
    d['Subject'] = TAG+' and VIC BC comparison metadata'
    d['Keywords'] = TAG+' VIC BC'
    d['CreationDate'] = datetime.datetime(2021, 4, 21)
    d['ModDate'] = datetime.datetime.today()