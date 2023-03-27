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
cd "/cluster/data6/menaka/AltiMaP"

#CaMA-Flood directory
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"

# map name
# map="glb_06min"
# map="amz_06min"
# map="glb_01min"
map="conus_06min"

glb_map="glb_06min"  # need to change according to map

# Higher resolution data
TAG="3sec"

# out put directory
outdir="./out"

# 
dataname="HydroWeb"

# obstxt
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"$glb_map"_20210909.txt"
obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"$glb_map"_20230327.txt"

# outtxt
day=$(date +"%Y%m%d")
# day='20210920'
outtxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"$map"_"$day".txt"

#
printf '%13s%64s%12s%12s%10s%17s%6s%12s%15s%10s%8s%8s%8s%14s%12s%12s%10s%8s%12s%10s\n' ID station dataname lon lat satellite flag elevation dist_to_mouth kx1 ky1 kx2 ky2 dist1 dist2 rivwth ix iy EGM08 EGM96 > tmp.txt
python src/regional_alloc.py $dataname $map $glb_map $CaMa_dir $obstxt $outtxt >> tmp.txt 
mv tmp.txt $outtxt
echo "$outtxt created." 