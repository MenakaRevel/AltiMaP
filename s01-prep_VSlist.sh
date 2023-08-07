#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=60gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N VS_list

#source ~/.bashrc

export OMP_NUM_THREADS=20

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/AltiMaP"

# `pwd`

# Data name
# dataname="HydroWeb"
dataname="Schneider2017"

# data directory
# datadir="/cluster/data6/menaka/HydroWeb/data"
datafile="./ESSD_paper/CryoSat2_Brahmaputra_list.txt"


# output directory
outdir="./inp"

python "./src/make_VSlist.py" $dataname $datafile $outdir