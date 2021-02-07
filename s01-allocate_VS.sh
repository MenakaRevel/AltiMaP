#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=100gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Alti_VS

#source ~/.bashrc

NCPUS=20
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
cd $PBS_O_WORKDIR

#CaMA-Flood directory
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"

# map name
map="glb_01min"

# Higher resolution data
TAG="3sec"

# out put directory
outdir="./out"

mkdir -p $outdir

USER=`whoami`

echo "            ID                                      station            dataname         lon       lat       ix      iy     ele_diff     EGM08     EGM96        satellite" > tmp.txt
SOUTH=-60
while [ $SOUTH -lt 90 ];
do
  WEST=-180
  while [ $WEST -lt 180 ];
  do
    CNAME=`./src/set_name $WEST $SOUTH`
    #echo $CNAME ${CaMa_dir}/map/${map}/${TAG}/${CNAME}.catmxy.bin
    if [ -f ${CaMa_dir}/map/${map}/${TAG}/${CNAME}.catmxy.bin ]; then
        for data in "HydroWeb" "CGLS" "HydroSat" "GRRATS"; # "ICESat";
        do
            flag=`python ./src/avalability_data.py $data $WEST $SOUTH`
            # echo $flag
            if [ $flag = 1 ]; then
                echo "./src/allocate_VS $WEST $SOUTH $data"
                ./src/allocate_VS $WEST $SOUTH $data $CaMa_dir $map $TAG $outdir >> tmp.txt &
                ## for parallel computation using multiple CPUs 
                NUM=`ps aux -U $USER | grep /src/allocate_VS | wc -l | awk '{print $1}'`
                #echo $USER $NUM
                while [ $NUM -gt $NCPUS ];
                do
                    sleep 1
                    NUM=`ps aux -U $USER | grep /src/allocate_VS | wc -l | awk '{print $1}'`
                    echo $USER $NUM
                done
            fi
        done
    fi
  WEST=$(( $WEST + 10 ))
  done
SOUTH=$(( $SOUTH + 10 ))
done


wait 
mv tmp.txt ./out/altimetry_${map}.txt