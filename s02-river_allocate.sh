#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=10gb
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

NCPUS=20
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Altimetry"

#CaMA-Flood directory
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"

# map name
# map="glb_01min"
map="glb_06min"

# Higher resolution data
TAG="3sec"

# date name
dataname="HydroWeb"

# out put directory
# outdir="/cluster/data6/menaka/Altimetry/results"
outdir="/cluster/data6/menaka/Altimetry/fig"

# observation list
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_test.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210531.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210602.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210618.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210701.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210706.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210709.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210714.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210808.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210817.txt"
# obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210909.txt"
obstxt="/cluster/data6/menaka/Altimetry/out/altimetry_"$map"_20210920.txt"

mkdir -p ${outdir}/${dataname}/high_res

USER=`whoami`

# python src/river_map.py $map $CaMa_dir $NCPUS & #> /dev/null 2>&1 & 
# python src/high_res.py "CONGO" "CONGO" $dataname $outdir $map $CaMa_dir $TAG $obstxt #&
# python src/high_res.py "AMAZONAS" "CONGO" $dataname $outdir $map $CaMa_dir $TAG $obstxt &
# python src/high_res.py "LENA" "LENA" $dataname $outdir $map $CaMa_dir $TAG $obstxt #&
#  "AMAZONAS" "CONGO"; #

# create figures for each river reach
python src/high_res.py "AMAZONAS" "MADEIRA" $dataname $outdir $map $CaMa_dir $TAG $obstxt &
python src/high_res.py "AMAZONAS" "PURUS" $dataname $outdir $map $CaMa_dir $TAG $obstxt &
python src/high_res.py "AMAZONAS" "NEGRO" $dataname $outdir $map $CaMa_dir $TAG $obstxt &
python src/high_res.py "AMAZONAS" "SOLIMOES" $dataname $outdir $map $CaMa_dir $TAG $obstxt &
python src/high_res.py "CONGO" "KASAI" $dataname $outdir $map $CaMa_dir $TAG $obstxt &
python src/high_res.py "GANGES-BRAHMAPUTRA" "GANGES" $dataname $outdir $map $CaMa_dir $TAG $obstxt &
python src/high_res.py "GANGES-BRAHMAPUTRA" "BRAHMAPUTRA" $dataname $outdir $map $CaMa_dir $TAG $obstxt &


for rivername in "AMAZONAS" "CONGO" "MEKONG" "VOLGA" "NIGER" "IRRAWADDY" "MISSISSIPPI" "LENA";  # "PARANA" "YANGTZE" "AMUR" "DANUBE" "ORINOCO" "PARNAIBA" "SAO-FRANCISCO"
do
    python src/high_res.py $rivername $rivername $dataname $outdir $map $CaMa_dir $TAG $obstxt &
    # echo "python src/high_res.py $rivername $rivername $dataname $outdir $map $CaMa_dir $TAG $obstxt" &
    ## for parallel computation using multiple CPUs 
    NUM=`ps aux -U $USER | grep src/high_res.py | wc -l | awk '{print $1}'`
    echo $USER $NUM
    while [ $NUM -gt $NCPUS ];
    do
        sleep 1
        NUM=`ps aux -U $USER | grep src/high_res.py | wc -l | awk '{print $1}'`
        echo $USER $NUM
    done
done
wait

conda deactivate