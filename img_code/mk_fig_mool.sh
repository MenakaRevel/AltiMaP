#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F40
#PBS -l select=1:ncpus=40:mem=180gb
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

NCPUS=40
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Altimetry/img_code"


mkdir -p fig

# python f04-flag_map_example.py & 

# python f06-unrealistic_map_dist.py &

python f08-comapre_expert_ordinary_map.py

wait

conda deactivate