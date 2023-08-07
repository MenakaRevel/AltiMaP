#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F10
#PBS -l select=1:ncpus=10:mem=180gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N mk_fig

#===========================
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

NCPUS=10
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/AltiMaP/ESSD_paper"


mkdir -p fig

python f04-flag_map_uparea.py & 

# python f05-biased_VS_map_dist.py &

# python f06-unrealistic_map_dist.py &

python f07-map_compare_methods.py &

wait

conda deactivate