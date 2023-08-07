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

#===========================
# import virtual environment
# source ~/.bashrc
# source ~/.bash_conda

# source activate pydef

which python

NCPUS=20
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/AltiMaP"

#CaMA-Flood directory
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"

# higher resolution tag
# TAG="15sec"
TAG="3sec"

# map name
map="glb_06min"
# map="amz_06min"

# date name
dataname="HydroWeb"
# dataname="CGLS"

# observation list
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"$map"_20221205.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"$map"_20230327.txt"
# obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"$map"_20230406.txt"
obstxt="/cluster/data6/menaka/AltiMaP/out/altimetry_"$map"_20230407.txt"

# out dir
outdir="./out"
mkdir -p $outdir

day=$(date +"%Y%m%d")
outname="biased_removed_altimetry_"$map"_"$day".txt"
# outname="biased_removed_altimetry_"$map"_20230407.txt"

# threshold for finding outliers
threshold=15.0  #10.0 #m

# method="static" # use threshold
method="dynamic" # use 3*std

printf '%13s%64s%12s%12s%10s%17s%6s%12s%15s%10s%8s%8s%8s%14s%12s%12s%10s%8s%12s%10s\n' ID station dataname lon lat satellite flag elevation dist_to_mouth kx1 ky1 kx2 ky2 dist1 dist2 rivwth ix iy EGM08 EGM96 > tmp.txt
python src/unreal_obs.py $dataname $map $CaMa_dir $TAG $obstxt $threshold $method >> tmp.txt  #& #> /dev/null 2>&1 & 

echo "saving output to $outdir/$outname"
mv "tmp.txt" $outdir/$outname

wait

# conda deactivate