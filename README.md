### Dataset Varaibles
# VS metadata
ID: Identification number of VS
station         :VS name
dataname        :dataset name
lon 	        :longitude [degrees east]
lat	            :latitude [degrees north]
satellite	    :name of the satellite
# MERIT Hydro-related	
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
# Coarse-resolution river network-related	
ix	            :x-coordinate with respect to coarse resolution
iy	            :y-coordinate with respect to coarse resolution
EGM08	        :EGM 2008 datum elevation[m]
EGM96	        :EGM 1996 datum elevation[m]


### Altimetry Virtual Stations on CaMa-Flood River Network

1. Allocate virtual stations to the CaMa-Flood map
```bash
sh s01-allocate_VS.sh
```
2. Check the river network allocation
```bash
s02-river_allocate.sh
```
3. Compare the common Virtual Station on the same CaMa-Flood grid
```bash
s03-commonVS_plot.sh 
``` 
4. Compare common Virtual Station in close proximity
```bash
s03-along_river_plot.sh
```

### Figure Codes
1. 