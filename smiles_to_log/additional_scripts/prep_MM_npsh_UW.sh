#!/bin/bash
export SCHRODINGER=/home/sheins/schrodinger2022-4

COM_REFERENCE=${PWD}'/mm_reference.com' #run a manual macromodel job and save the .com file after removing the first two lines with input/output file names

if [[ ! -d ${PWD}/output/ ]]; then #if the output directory doesn't exist, make it y
mkdir ${PWD}/output/
fi

if [[ ! -d ${PWD}/output/clustered ]]; then #if the output directory doesn't have a cluster folder, make one
mkdir ${PWD}/output/clustered
fi

for file in $PWD/*.sdf; do #convert all .sdfs to maestro files
echo $file
filename=$(basename $file) #file.sdf
nopathnoext=${filename%.*} #file

echo "convert to .mae"
echo " "
$SCHRODINGER/utilities/structconvert ${nopathnoext}.sdf ${nopathnoext}.mae

mkdir $nopathnoext #make a directory for each compound for all output files and move the maestro file into it
mv ${nopathnoext}.mae ${nopathnoext}/${nopathnoext}.mae

cd ${nopathnoext}/

echo "${nopathnoext}.mae" > ${nopathnoext}.com #make the file names consistent on all lines
echo "${nopathnoext}_out.maegz" >> ${nopathnoext}.com
cat $COM_REFERENCE >> ${nopathnoext}.com #for each maestro file, add names for input and output into the reference com file, ultimately generating a com file for each compound

$SCHRODINGER/bmin -WAIT -HOST "localhost:1" ${nopathnoext}

cd ../

done
