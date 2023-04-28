### Dataset Variables
#### VS metadata
ID: Identification number of VS
station         :VS name
dataname        :dataset name
lon 	        :longitude [degrees east]
lat	            :latitude [degrees north]
satellite	    :name of the satellite
#### MERIT Hydro-related	
flag	        :allocation flag
elevation	    :elevation at VS location on MERIT Hydro [m]
dist_to_mouth	:distance to the unit-catchment mouth [km]
kx1	            :best x-coordinate with respect to the 10°×10° higher resolution tile
ky1	            :best y-coordinate with respect to the 10°×10° higher resolution tile
kx2	            :second-best option of x-coordinate with respect to the 10°×10° high-resolution tile
ky2	            :second-best option of y-coordinate with respect to the 10°×10° high-resolution tile
dist1	        :distance from the second-best location to the VS [km]
dist2	        :distance from the second-best location to the VS [km]
rivwth	        :River width of the allocated location [m]
#### Coarse-resolution river network-related	
ix	            :x-coordinate with respect to coarse resolution
iy	            :y-coordinate with respect to coarse resolution
EGM08	        :EGM 2008 datum elevation[m]
EGM96	        :EGM 1996 datum elevation[m]


### Mapping Altimetry Virtual Stations on CaMa-Flood River Network

1. Prepare virtual stations list with geoid conversation
```bash
sh s01-prep_VSlist.sh
```

2. Allocate virtual stations to the CaMa-Flood map
```bash
sh s02-allocate_VS.sh
```
4. Identify largely biased virtual stations
```bash
sh s03-unreal_obs.sh
```

### Input data
1. list of virtual station
2. MERIT-Hydro river parameters
3. CaMa-Flood river network map

#### Input data style for s01-prep_VSlist.sh
filename={datadir}/{dataname}_VS
identifier,river,basin,country,satellite,track_nb,start_date,end_date,latitude,longitude,status,validation

#### Input data style for s02-allocate_VS.sh
filename={dataname}_Station_list.txt
ID | Station | River | Basin | Country | lon | lat | elevation | EGM08 | EGM96 | Satellite | Start Date | End Date | Status

#### Input data style for s03-unreal_obs.sh
output from s02-allocate_VS.sh