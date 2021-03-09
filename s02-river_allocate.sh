#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F10
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Riv_Map

#source ~/.bashrc

# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

NCPUS=10
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
cd $PBS_O_WORKDIR

#CaMA-Flood directory
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"

# map name
map="glb_01min"
# map="glb_06min"

# Higher resolution data
TAG="3sec"

# out put directory
outdir="./fig/river_network"

mkdir -p $outdir
mkdir -p tmp

# USER=`whoami`

python src/river_map.py $map $CaMa_dir $NCPUS & #> /dev/null 2>&1 & 

wait

rm -r tmp
conda deactivate