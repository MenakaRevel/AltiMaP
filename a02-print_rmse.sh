#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F10
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N RMSE_VS

#source ~/.bashrc

NCPUS=10
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Altimetry"

#CaMA-Flood directory
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"

# map name
map="glb_06min"
# map="amz_06min"
# map="glb_01min"

glb_map="glb_06min"  # need to change according to map

# Higher resolution data
TAG="3sec"

# out put directory
outdir="./out"

# 
dataname="HydroWeb"

# obstxt
obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$glb_map"_20210909.txt"

# NetCDF file
netcdf='/cluster/data6/menaka/Altimetry/results/HydroWeb/hydroweb-hydroda_cmf_daily_wse_VIC_BC.nc'

printf '%67s%10s%10s\n' Station RMSE BIAS > tmp.txt

python src/rmse.py >> tmp.txt

# outtxt
day=$(date +"%Y%m%d")
outtxt="/cluster/data6/menaka/Altimetry/out/rmse_VIC_BC_"$map"_"$day".txt"
mv tmp.txt $outtxt