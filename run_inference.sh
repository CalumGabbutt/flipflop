#!/bin/sh

datafile=data/beta_ticktock.csv
patientinfofile=data/patientinfo.csv
outputdir=examples
samplename=PA4
Smin=3
Smax=7
nlive=300

mkdir -p ${outputdir}

for ((S=$Smin;i<=$Smax;i++)); 
do
    python3 inference.py $datafile $patientinfofile $outputdir $samplename $S --nlive $nlive --verbose
done 
python3 combine_samples.py $datafile $patientinfofile $outputdir $samplename
