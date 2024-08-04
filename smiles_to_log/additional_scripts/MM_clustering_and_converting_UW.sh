#!/bin/bash

export SCHRODINGER=/home/sheins/schrodinger2022-4
WORKDIR=$PWD

if [[ ! -d ${PWD}/output/ ]]; then #usually already made from the prep_MM_npsh.bash
mkdir ${PWD}/output/
fi

if [[ ! -d ${PWD}/output/clustered ]]; then #if the output directory doesn't have a cluster folder, make one
mkdir ${PWD}/output/clustered
fi

for PDB in $WORKDIR/*/;do #for every folder in working directory

cd $PDB #enter the folder

for file in $PDB/*.log;do #for every .log file in that folder
filename=$(basename $file) #save the file name as the variable we'll use
nopathnoext=${filename%.*}

if [[ -f ${nopathnoext}.log ]];then #if the log file exists, continue
lines=`grep -A 1 "Final report; processing .tmp file:" ${nopathnoext}.log` #look for the final report in the log file, and save it to the lines variable
echo ${nopathnoext} has ${lines//[^0-9]/} conformers #use that information to print how many conformers that compound has
number=${lines//[^0-9]/} #store the number of conformers (the number in that line)

if [ $number -le 20 ] #if there <= 20 total conformers, move/rename the conformational search file output to the acid output folder
then
echo 'converting '${nopathnoext}'.maegz to .sdf'
$SCHRODINGER/utilities/structconvert ${nopathnoext}_out.maegz ${nopathnoext}_conf.sdf
mv ${nopathnoext}_conf.sdf ../output/${nopathnoext}_conf.sdf

else #if there are more than 20 conformers, move/rename to the clustering folder
$SCHRODINGER/run conformer_cluster.py -j ${nopathnoext}_cluster -n 0  -r -rep -WAIT ${nopathnoext}_out.maegz #if you want a certain number of clusters, use "-n 5", 0 uses the Kelley Penalty minimum

echo 'converting '${nopathnoext}'.maegz to .sdf'
$SCHRODINGER/utilities/structconvert ${nopathnoext}_cluster_ligand1_representatives.maegz ${nopathnoext}_clust.sdf
$SCHRODINGER/utilities/structconvert ${nopathnoext}_out.maegz ${nopathnoext}_conf.sdf
mv ${nopathnoext}_clust.sdf ../output/${nopathnoext}_clust.sdf
mv ${nopathnoext}_conf.sdf ../output/clustered/${nopathnoext}_conf.sdf

fi

fi
done

cd $WORKDIR

done

