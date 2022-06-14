#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F10
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Reg_VS

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
# map="glb_06min"
map="amz_06min"
# map="glb_01min"

glb_map="glb_06min"  # need to change according to map

# Higher resolution data
TAG="3sec"

# out put directory
outdir="./out"

# 
dataname="HydroWeb"

# obstxt
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$glb_map"_20210909.txt"
obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$glb_map"_20210920.txt"

# outtxt
day=$(date +"%Y%m%d")
# day='20210920'
outtxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_"$day".txt"

#
python src/regional_alloc.py $dataname $map $glb_map $CaMa_dir $obstxt $outtxt &