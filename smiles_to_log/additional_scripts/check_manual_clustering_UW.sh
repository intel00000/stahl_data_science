#!/bin/bash

WORKDIR=$PWD

for i in *.sdf; do babel -isdf $i -ogau ${i%.*}-.com -m; done

for file in *.sdf;do #for every .sdf file in that folder
filename=$(basename $file) #save the file name as the variable we'll use
nopathnoext=${filename%.*}

if [[ $filename == *"clust"* ]];then
if [[ ! -f ${nopathnoext}-2.com ]];then
echo $filename
rm ${nopathnoext}.sdf
rm ${nopathnoext}-1.com
fi
fi

done
