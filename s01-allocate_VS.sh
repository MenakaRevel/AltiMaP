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
# cd $PBS_O_WORKDIR
cd "/cluster/data6/menaka/Altimetry"

#CaMA-Flood directory
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"

# map name
map="glb_06min"
# map="amz_06min"
# map="glb_01min"

glb_map="glb_06min"  # need to change according to map

# Higher resolution data
TAG="3sec"

# out put directory
outdir="./out"

mkdir -p $outdir

USER=`whoami`

echo "making regional bifurication tags"
fbiftag="${CaMa_dir}/map/${map}/biftag.bin"
if [ ! -f ${fbiftag} ]; then
  ./src/regional_biftag $map $glb_map $CaMa_dir
fi

echo "starting calculations........"
# echo "            ID                                      station            dataname         lon       lat       ix      iy     ele_diff     EGM08     EGM96        satellite" > tmp.txt
printf '%30s%67s%12s%12s%10s%10s%8s%12s%10s%10s%17s%15s%6s%10s%8s%8s%8s%12s%12s\n' ID station dataname lon lat ix iy elevation EGM08 EGM96 satellite dist_to_mouth flag kx1 ky1 kx2 ky2 dist1 dist2 > tmp.txt
SOUTH=-60
while [ $SOUTH -lt 90 ];
do
  WEST=-180
  while [ $WEST -lt 180 ];
  do
    CNAME=`./src/set_name $WEST $SOUTH`
    # echo $CNAME #${CaMa_dir}/map/${map}/${TAG}/${CNAME}.catmxy.bin
    if [ -s ${CaMa_dir}"/map/"${map}/${TAG}/${CNAME}".catmxy.bin" ]; then
        for data in "HydroWeb"; #"CGLS" "HydroSat" "GRRATS"; # "ICESat";
        do
            flag=`python ./src/avalability_data.py $data $WEST $SOUTH`
            # echo $flag
            if [ $flag = 1 ]; then
                echo "./src/allocate_VS $WEST $SOUTH $data"
                ./src/allocate_VS $WEST $SOUTH $data $CaMa_dir $map $TAG $outdir >> tmp.txt &
                ## for parallel computation using multiple CPUs 
                NUM=`ps aux -U $USER | grep /src/allocate_VS | wc -l | awk '{print $1}'`
                # echo $USER $NUM
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

# day=`printf '%(%Y%m%d)T\n' -1`
day=$(date +"%Y%m%d")
mv tmp.txt ${outdir}/altimetry_${map}_${day}.txt
echo "Saving ..."
echo ${outdir}/altimetry_${map}_${day}.txt