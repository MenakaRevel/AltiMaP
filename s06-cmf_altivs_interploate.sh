#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=100gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N MK_netCDF_inter

#source ~/.bashrc
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

NCPUS=20
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Altimetry"

#CaMA-Flood directory
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"

# current working directory
curr_pwd='/cluster/data6/menaka/Altimetry'

#Simulation period
syear=2002
eyear=2013

# runoff_folder = '/cluster/data6/menaka/ensemble_org/CaMa_out/GLBVIC_BC_USED/'
runoff_folder='/cluster/data6/menaka/ensemble_org/CaMa_out/GLBVIC_BC001'

#observation dataset
TAG='HydroWeb'

out_pwd=${curr_pwd}"/results/"${TAG}

mkdir -p ${out_pwd}

python src/CMF_wse_inter.py $TAG $syear $eyear $runoff_folder $CaMa_dir $curr_pwd $out_pwd &

wait

conda deactivate