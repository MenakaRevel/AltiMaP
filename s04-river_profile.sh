#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=20gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Unreal_VS

# # import virtual environment
# source ~/.bashrc
# source ~/.bash_conda

# source activate pydef

# which python

NCPUS=20
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Altimetry"

# out dir
outdir="./out"
mkdir -p $outdir

#CaMA-Flood directory
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"

# high resolution tag
# TAG="15sec"
TAG="3sec"

# map name
map="glb_06min"

# date name
dataname="HydroWeb"

# observation list
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_test.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210531.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210602.txt"
obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210618.txt"

# unreal obseravation list
notxt="/cluster/data6/menaka/Altimetry/out/unreal_obs_20210622.txt"

# # out put figure directory
# figdir="./fig"

# mkdir -p $figdir
day=$(date +"%Y%m%d")

python src/river_profile.py $map $TAG $notxt $obstxt $CaMa_dir  > "$outdir/intporation_obs_${day}.txt" #& #> /dev/null 2>&1 & 

# wait

# conda deactivate