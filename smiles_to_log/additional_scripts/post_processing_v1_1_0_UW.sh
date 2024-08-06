#!/bin/bash

WORKDIR=$PWD

rm *.chk
rm *.fchk
rm slurm*
rm fort.7
find . -type d -empty -delete

if [[ ! -d ${PWD}/resubmit/ ]]; then 
mkdir ${PWD}/resubmit/
fi

if [[ ! -d ${PWD}/logs/ ]]; then
mkdir ${PWD}/logs/
fi

if [[ ! -d ${PWD}/term_error/ ]]; then
mkdir ${PWD}/term_error/
fi

if [[ ! -d ${PWD}/imag/ ]]; then
mkdir ${PWD}/imag/
fi

for file in *.com;do #for every .com file in that folder
filename=$(basename $file) #save the file name as the variable we'll use
nopathnoext=${filename%.*}

if [[ -f ${nopathnoext}.log ]];then #if the log file exists, move it to the logs output
mv ${nopathnoext}.log logs/${nopathnoext}.log
else

mv ${nopathnoext}.com resubmit/${nopathnoext}.com
fi

done

cd logs/
python /home/sheins/log_check_for_processing_v1_1_0.py
