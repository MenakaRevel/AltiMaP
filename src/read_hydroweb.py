#!/usr/bin/python
import os
import re
import datetime
import numpy as np
import pandas as pd ## convert the function to pandas @2023/-7
################################
# Fuctions to read HydroWeb data
################################
def hydroweb_river_name(mapname="glb_15min"):
    # directory
    fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    df=pd.read_csv(fname, sep='\s+', header=0)
    df=df.dropna()
    #--
    # Split the station column
    df['river'] = df['station'].apply(lambda x: re.split("_", x)[1])

    # Convert the 'river' column to a list
    river = df['river'].tolist()
    return river
################################
def get_hydroweb(mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    # fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    df=pd.read_csv(fname, sep='\s+', header=0)
    df=df.dropna()
    #--
    # Split the station column and create new columns
    df['river'] = df['station'].apply(lambda x: re.split("_", x)[1])
    df['pname'] = df['station']
    df['xlist'] = df['ix'].apply(lambda x: int(x) - 1)
    df['ylist'] = df['iy'].apply(lambda x: int(x) - 1)
    df['egm_d'] = df['EGM96'].astype(float) - df['EGM08'].astype(float)

    # Convert the new columns to lists
    river = df['river'].tolist()
    pname = df['pname'].tolist()
    xlist = df['xlist'].tolist()
    ylist = df['ylist'].tolist()
    egm_d = df['egm_d'].tolist()
    return river,pname,xlist,ylist,egm_d
##################################
def get_hydroweb_loc(rivername,mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    df=pd.read_csv(fname, sep='\s+', header=0)
    df=df.dropna()
    #--
    # Split the station column and create new columns
    df['river'] = df['station'].apply(lambda x: re.split("_", x)[1])
    df['pname'] = df['station']
    df['xlist'] = df['ix'].apply(lambda x: int(x) - 1)
    df['ylist'] = df['iy'].apply(lambda x: int(x) - 1)
    df['EGM08'] = df['EGM08'].astype(float)
    df['EGM96'] = df['EGM96'].astype(float) 

    # Filter rows where rivername equals riv
    df_filtered = df[df['river'] == rivername]

    # Convert the new columns to lists
    pname = df_filtered['pname'].tolist()
    xlist = df_filtered['xlist'].tolist()
    ylist = df_filtered['ylist'].tolist()
    egm08 = df_filtered['EGM08'].tolist()
    egm96 = df_filtered['EGM96'].tolist()

    return pname,xlist,ylist,egm08,egm96
##################################
def get_hydroweb_locs(mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    # directory
    df=pd.read_csv(fname, sep='\s+', header=0)
    df=df.dropna()
    # Split the station column and create new columns
    df['ID']            = df['ID']
    df['river']         = df['station'].apply(lambda x: re.split("_", x)[1])
    df['station']       = df['station']
    df['lon']           = df['lon'].astype(float)
    df['lat']           = df['lat'].astype(float)
    df['ix']            = df['ix'].apply(lambda x: int(x) - 1)
    df['iy']            = df['iy'].apply(lambda x: int(x) - 1)
    df['elevation']     = df['elevation'].astype(float)
    df['EGM08']         = df['EGM08'].astype(float)
    df['EGM96']         = df['EGM96'].astype(float)
    df['satellite']     = df['satellite'].str.strip()
    df['dist_to_mouth'] = df['dist_to_mouth'].astype(float)
    df['flag']          = df['flag'].astype(int)

    # Convert the new columns to lists
    nums  = df['ID'].tolist()
    river = df['river'].tolist()
    pname = df['station'].tolist()
    lons  = df['lon'].tolist()
    lats  = df['lat'].tolist()
    xlist = df['ix'].tolist()
    ylist = df['iy'].tolist()
    elev  = df['elevation'].tolist()
    egm08 = df['EGM08'].tolist()
    egm96 = df['EGM96'].tolist()
    llsat = df['satellite'].tolist()
    ldtom = df['dist_to_mouth'].tolist()
    lflag = df['flag'].tolist()
    return nums,river,pname,lons,lats,xlist,ylist,egm08,egm96,llsat,elev,ldtom,lflag
################################
def HydroWeb_WSE(station,syear,eyear,smon=1,emon=12,sday=1,eday=31,egm08=0.0,egm96=0.0,odir0="/cluster/data6/menaka/HydroWeb/data"):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)

    # Read data from the file using pandas
    fname=dir0+"/hydroprd_"+station+".txt"

    # Assuming the file has columns: 'Date', 'Time', 'WSE', 'Uncertainty'
    df = pd.read_csv(fname, skiprows=33, sep='\s+', names=['Date', 'Time', 'WSE', 'Uncertainty'])

    # Combine 'Date' and 'Time' columns into a single 'DateTime' column
    df['DateTime'] = pd.to_datetime(df['Date'] + ' ' + df['Time'])

    # Calculate the water surface elevation with the given adjustments
    df['WSE'] = df['WSE'] + egm08 - egm96

    # Filter the data within the specified date range
    df = df[(df['DateTime'] >= start) & (df['DateTime'] <= end)]

    # Calculate the time in days relative to the 'start' date
    df['Days'] = (df['DateTime'] - start).dt.days

    # Drop unnecessary columns
    df = df.drop(columns=['Date', 'Time', 'Uncertainty'])

    data = df['WSE'].values
    time = df['Days'].values

    return time, data
#####################################
def HydroWeb_continous_WSE(station,syear=2002,smon=10,sday=1,eyear=2020,emon=12,eday=31,egm08=0.0,egm96=0.0,dir0="/cluster/data6/menaka/HydroWeb/data"):
    #
    start=datetime.date(syear,smon,sday)
    end=datetime.date(eyear,emon,eday)

    # Read data from the file using pandas
    fname=dir0+"/hydroprd_"+station+".txt"

    # Assuming the file has columns: 'Date', 'Time', 'WSE', 'Uncertainty'
    df = pd.read_csv(fname, skiprows=33, sep='\s+', names=['Date', 'Time', 'WSE', 'Uncertainty'])

    # Combine 'Date' and 'Time' columns into a single 'DateTime' column
    df['DateTime'] = pd.to_datetime(df['Date'] + ' ' + df['Time'])

    # Calculate the water surface elevation with the given adjustments
    df['WSE'] = df['WSE'] + egm08 - egm96

    # Filter the data within the specified date range
    df = df[(df['DateTime'] >= start) & (df['DateTime'] <= end)]

    # Calculate the time in days relative to the 'start' date
    df['Days'] = (df['DateTime'] - start).dt.days

    # Create an array with WSE values, initialized with -9999.0 for no observations
    time = int((end - start).days) + 1
    data = np.ones(time, dtype=np.float32) * -9999.0

    # Calculate the corresponding indices for the 'data' array using vectorized operations
    time_indices = (df['Days']).values

    # Calculate the adjusted WSE values in bulk using numpy array operations
    wse_values = df['WSE'].values
    data[time_indices] = wse_values + egm08 - egm96

    return data
#####################################
def altimetry(name,mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
    # directory
    hydroweb="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
    #--
    f=open(fname,"r")
    lines=f.readlines()
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        station = line[1]
        if station==name:
            alti= float(line[6])
    return alti
################################
def metadata():
    # directory 
    hydroweb="/cluster/data6/menaka/HydroWeb/HydroWebStation_list.txt"
    #--
    f=open(hydroweb,"r")
    lines=f.readlines()
    sta={}
    for line in lines[1::]:
        line    = filter(None,re.split(" ",line))
        id      = int(line[0])
        station = line[1]
        river   = line[2]
        basin   = line[3]
        country = line[4]
        sat     = line[10]
        startdt = line[11]
        enddt   = line[13]
        try:
            status  = line[15].split("\n")[0]
        except:
            status = "N/A"
        #----
        sta[station]=[river,basin,country,sat,startdt,enddt,status]
    #----
    return sta
################################










#=========================================
## old codes
#=========================================
# # #!/usr/bin/python
# # import os
# # import re
# # import datetime
# # import numpy as np
# # import pandas as pd ## convert the function to pandas @2023/-7
# # ################################
# # # Fuctions to read HydroWeb data
# # ################################
# # def hydroweb_river_name(mapname="glb_15min"):
# #     # directory
# #     hydroweb="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
# #     #--
# #     river=[]
# #     #--
# #     f=open(hydroweb,"r")
# #     lines=f.readlines()
# #     for line in lines[1::]:
# #         line    = filter(None,re.split(" ",line))
# #         print line
# #         station = line[1]
# #         riv     = re.split("_",station)[1]
# #         river.append(riv)
# #     return river
# # ################################
# # def get_hydroweb(mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
# #     #--
# #     river=[]
# #     pname=[]
# #     xlist=[]
# #     ylist=[]
# #     egm_d=[]
# #     #--
# #     f=open(fname,"r")
# #     lines=f.readlines()
# #     for line in lines[1::]:
# #         line    = filter(None,re.split(" ",line))
# #         print line
# #         station = line[1]
# #         riv     = re.split("_",station)[1]
# #         ix      = int(line[4])-1
# #         iy      = int(line[5])-1
# #         EGM08   = float(line[8])
# #         EGM96   = float(line[9])
# #         river.append(riv)
# #         pname.append(station)
# #         xlist.append(ix)
# #         ylist.append(iy)
# #         egm_d.append(EGM96-EGM08)
# #     return river,pname,xlist,ylist,egm_d
# # ##################################
# # def get_hydroweb_loc(rivername,mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
# #     # directory
# #     # fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
# #     #--
# #     #river=[]
# #     pname=[]
# #     xlist=[]
# #     ylist=[]
# #     egm08=[]
# #     egm96=[]
# #     #--
# #     f=open(fname,"r")
# #     lines=f.readlines()
# #     for line in lines[1::]:
# #         line    = filter(None,re.split(" ",line))
# #         #print line
# #         station = line[1]
# #         riv     = re.split("_",station)[1]
# #         ix      = int(line[4])-1
# #         iy      = int(line[5])-1
# #         EGM08   = float(line[8])
# #         EGM96   = float(line[9])
# #         if rivername==riv:
# #             #river.append(riv)
# #             pname.append(station)
# #             xlist.append(ix)
# #             ylist.append(iy)
# #             egm08.append(EGM08)
# #             egm96.append(EGM96)
# #     return pname,xlist,ylist,egm08,egm96
# # ##################################
# # def get_hydroweb_locs(mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
# #     # directory
# #     # fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
# #     # fname="./out/altimetry_"+mapname+"_test.txt"
# #     #--
# #     nums=[]
# #     river=[]
# #     pname=[]
# #     lons =[]
# #     lats =[]
# #     xlist=[]
# #     ylist=[]
# #     leled=[]
# #     egm08=[]
# #     egm96=[]
# #     llsat=[]
# #     ldtom=[]
# #     lflag=[]
# #     #--
# #     f=open(fname,"r")
# #     lines=f.readlines()
# #     for line in lines[1::]:
# #         line    = filter(None,re.split(" ",line))
# #         #print line
# #         num     = line[0]
# #         station = line[1]
# #         riv     = re.split("_",station)[1]
# #         lon     = float(line[3])
# #         lat     = float(line[4])
# #         ix      = int(line[5])-1
# #         iy      = int(line[6])-1
# #         eled    = float(line[7])
# #         EGM08   = float(line[8])
# #         EGM96   = float(line[9])
# #         sat     = line[10].strip()
# #         dist    = float(line[11])
# #         flag    = int(line[12])
# #         #-----------------------
# #         nums.append(num)
# #         river.append(riv)
# #         pname.append(station)
# #         lons.append(lon)
# #         lats.append(lat)
# #         xlist.append(ix)
# #         ylist.append(iy)
# #         leled.append(eled)
# #         egm08.append(EGM08)
# #         egm96.append(EGM96)
# #         llsat.append(sat)
# #         ldtom.append(dist)
# #         lflag.append(flag)
# #     return nums,river,pname,lons,lats,xlist,ylist,egm08,egm96,llsat,leled,ldtom,lflag
# # ################################
# # def HydroWeb_WSE(station,syear,eyear,smon=1,emon=12,sday=1,eday=31,egm08=0.0,egm96=0.0,fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
# #     #
# #     start=datetime.date(syear,smon,sday)
# #     end=datetime.date(eyear,emon,eday)
# #     # read hydroweb
# #     #station="R_con_con_env_0429_01"
# #     #satellite=station.split("_")[2]
# #     #fname="/home/yamadai/data/Altimetry/HydroWeb_LEGOS/River/R_"+station
# #     # fname="/cluster/data6/menaka/HydroWeb/data/hydroprd_"+station+".txt"
# #     f=open(fname,"r")
# #     lines=f.readlines()
# #     f.close()
# #     head=33
# #     #--
# #     time=[] # time in days
# #     data=[] # WSE in [m]
# #     for line in lines[head::]:
# #         if line[0][0] == "#":
# #             continue
# #         line = re.split(" ",line)
# #         date = line[0]
# #         date = re.split("-",date)
# #         yyyy = int(date[0])
# #         mm   = int(date[1])
# #         dd   = int(date[2])
# #         if float(line[2]) >= 9999.0:
# #             continue
# #         wse  = float(line[2])+egm08-egm96
# #         #print yyyy, mm, dd, wse
# #         now  = datetime.date(yyyy,mm,dd)
# #         if now < start and now > end:
# #             continue
# #         data.append(wse)
# #         lag  = int((now-start).days)
# #         time.append(lag)
# #     return time, data
# # #####################################
# # def HydroWeb_continous_WSE(station,syear=2002,smon=10,sday=1,eyear=2020,emon=12,eday=31,egm08=0.0,egm96=0.0,dir0="/cluster/data6/menaka/HydroWeb/data"):
# #     #
# #     start=datetime.date(syear,smon,sday)
# #     end=datetime.date(eyear,emon,eday)
# #     # read hydroweb
# #     #station="R_con_con_env_0429_01"
# #     #satellite=station.split("_")[2]
# #     #fname="/home/yamadai/data/Altimetry/HydroWeb_LEGOS/River/R_"+station
# #     # fname="/cluster/data6/menaka/HydroWeb/data/"hydroprd_"+station+".txt"
# #     fname=dir0+"/hydroprd_"+station+".txt"
# #     f=open(fname,"r")
# #     lines=f.readlines()
# #     f.close()
# #     head=33
# #     #--
# #     time=int((end-start).days + 1) # time in days
# #     data=np.ones([time],np.float32)*-9999.0 # WSE in [m] # -9999.0 for no observations
# #     for line in lines[head::]:
# #         if line[0][0] == "#":
# #             continue
# #         line = re.split(" ",line)
# #         date = line[0]
# #         date = re.split("-",date)
# #         yyyy = int(date[0])
# #         mm   = int(date[1])
# #         dd   = int(date[2])
# #         if float(line[2]) >= 9999.0:
# #             continue
# #         wse  = float(line[2])
# #         #print yyyy, mm, dd, wse
# #         now  = datetime.date(yyyy,mm,dd)
# #         #print yyyy, mm, dd, int((now-start).days), wse
# #         if now > start and now < end:
# #             nowtime=int((now-start).days)
# #             data[nowtime]=wse+egm08-egm96
# #     return data
# # #####################################
# # def altimetry(name,mapname="glb_15min",fname="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_glb_15min.txt"):
# #     # directory
# #     hydroweb="/cluster/data6/menaka/HydroWeb/HydroWeb_alloc_"+mapname+".txt"
# #     #--
# #     f=open(fname,"r")
# #     lines=f.readlines()
# #     for line in lines[1::]:
# #         line    = filter(None,re.split(" ",line))
# #         station = line[1]
# #         if station==name:
# #             alti= float(line[6])
# #     return alti
# # ################################
# # def metadata():
# #     # directory 
# #     hydroweb="/cluster/data6/menaka/HydroWeb/HydroWebStation_list.txt"
# #     #--
# #     f=open(hydroweb,"r")
# #     lines=f.readlines()
# #     sta={}
# #     for line in lines[1::]:
# #         line    = filter(None,re.split(" ",line))
# #         id      = int(line[0])
# #         station = line[1]
# #         river   = line[2]
# #         basin   = line[3]
# #         country = line[4]
# #         sat     = line[10]
# #         startdt = line[11]
# #         enddt   = line[13]
# #         try:
# #             status  = line[15].split("\n")[0]
# #         except:
# #             status = "N/A"
# #         #----
# #         sta[station]=[river,basin,country,sat,startdt,enddt,status]
# #     #----
# #     return sta
# # ################################