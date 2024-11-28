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
# dataname="Schneider2017"
# dataname="Prakatgauge"
dataname="SWOTMackenzie"

# data directory
# datadir="/cluster/data6/menaka/HydroWeb/data"
# datadir="/cluster/data6/menaka/HydroWeb/data_v2023"
# datadir="/work/a06/menaka/Prakat_Model_Scale/inp"
datadir="/cluster/data6/menaka/AltiMaP/inp"

# datafile
# datafile="./ESSD_paper/CryoSat2_Brahmaputra_list.txt"
# datafile="list_hydroprd_202309151209_rivers.csv"
# datafile="gauge_Amazon_list.txt"
datafile="SWOT_Mackenzie_Station_list.txt"


# output directory
outdir="./inp"

python "./src/make_VSlist.py" $dataname $datafile $datadir $outdir